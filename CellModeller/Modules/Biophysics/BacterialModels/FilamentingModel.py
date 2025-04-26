from typing import Iterable, Optional
from CellModeller.Biophysics.BacterialModels.Bacterium2D import Bacterium2D

import numpy as np

from CellModeller.CellState import CellState


class BacteriumFilamenting(Bacterium2D):
    @staticmethod
    def renderer_attrs(
        max_cells: int,
        max_contacts: int,
        alignment: int,
    ) -> dict[str, tuple[str, np.dtype]]:
        del max_contacts, alignment
        return {
            "RodCells": (
                "BacteriumFilamenting",
                np.dtype(
                    {
                        "names": [
                            "position",
                            "direction",
                            "radii+cellLength+active",
                            "cellColor",
                        ],
                        "formats": ["4f", "4f", "4f", "4f"],
                        "offsets": [
                            0,
                            max_cells * 4 * 4,
                            2 * max_cells * 4 * 4,
                            3 * max_cells * 4 * 4,
                        ],
                    }
                ),
            )
        }

    def __init__(
        self,
        jitter_z: bool = False,
        bin_width: int = 10,
        cgs_tol: float = 0.000001,
        time_factor: float = 400,
        max_iter: Optional[int] = None,
        conlim: Optional[float] = None,
    ):
        super().__init__(jitter_z, bin_width, cgs_tol, time_factor, max_iter, conlim)
        add_numpy_attrs = np.dtype(
            [
                ("filamenting", "?"),
                ("joint", "?"),
                ("filaments", "2i"),
                ("filamentNeighbours", "2i"),
                ("joint_growth", "f"),
            ]
        )
        if self.numpy_attrs is not None:
            self.numpy_attrs = np.dtype(self.numpy_attrs.descr + add_numpy_attrs.descr)

    def buildMatrices(self, cell_states: Iterable[CellState]) -> None:
        neighbourState = next(iter(cell_states))
        for cellState in cell_states:
            neighbourhood = tuple(cellState.neighbourhood)

            if cellState.joint:
                f1index, f2index = cellState.filaments
                f1 = self.neighbourhoods[neighbourhood][f1index]
                f2 = self.neighbourhoods[neighbourhood][f2index]

                contact1, overlap_dist1, contact2, overlap_dist2 = (
                    self.find_contact_filament(cellState, f1, f2)
                )
                # print(contact1, "\n", contact2, "\n", overlap_dist1, overlap_dist2)

                self.build_B_filament(cellState, f1, contact1, overlap_dist1)
                self.build_B_filament(cellState, f2, contact2, overlap_dist2)

                # print(cellState.time)
                self.build_B_turn(f1, f2, cellState.pos)
                # print(cellState.time, contact, spring_force)

                # effective midpoint of filaments for weight calculations
                # can multiply by 1e9 to fix joint
                cellState.length = (f1.length + f2.length) / 2
                self.build_M(cellState)
                cellState.length = 0
            else:
                self.build_M(cellState)

            for neighbourIndex in self.neighbourhoods[neighbourhood]:
                if cellState.index < neighbourIndex:
                    neighbourState.change_view(neighbourIndex)
                    if not (
                        cellState.filamenting
                        and (
                            neighbourState.index in cellState.filaments
                            or neighbourState.index in cellState.filamentNeighbours
                        )
                    ):
                        contact, overlap_dist = self._find_contact(
                            cellState, neighbourState
                        )
                        if contact is not None:
                            self.build_B(
                                cellState, neighbourState, contact, overlap_dist
                            )

    def build_B_turn(
        self,
        f1: CellState,
        f2: CellState,
        center: np.ndarray,
    ) -> None:
        a1 = f1.pos - center
        a2 = f2.pos - center
        a1 = a1 / np.linalg.norm(a1)
        a2 = a2 / np.linalg.norm(a2)
        # n1 = n1 / np.linalg.norm(n1)
        # n2 = n2 / np.linalg.norm(n2)
        n1 = f1.dir * np.sign(np.dot(f1.dir, a1))
        n2 = f2.dir * np.sign(np.dot(f2.dir, a2))
        spring_force = np.linalg.norm(
            np.cross(n1, n2) / (1 - np.dot(n1, n2)) * (1 + np.dot(n1, n2))
        )
        overlap_dist = spring_force
        # FIX: why is this unstable?!??!!
        # print(n1, n2, spring_force)

        i = f1.index
        j = f2.index
        self.B_build[0] += (
            list(-a1 - np.cross(n1, np.cross(n1, n2)))
            + list(np.cross(np.cross(n1, n2), f1.dir))
            # + [np.dot(contact[0], joint.dir) * np.dot(contact[2], joint.dir)]
            + [0]
        )
        self.B_build[1] += [self.contact_count] * 7
        self.B_build[2] += list(range(7 * i, 7 * i + 7))

        self.B_build[0] += (
            list(-a2 - np.cross(n2, np.cross(n2, n1)))
            + list(
                np.cross(
                    np.cross(n2, n1),
                    f2.dir,
                )
                # + np.cross(filament.dir, np.cross(other.dir, filament.dir))
                # / filament.length
            )
            # + [np.dot(-contact[1], filament.dir) * np.dot(contact[3], filament.dir)]
            + [0]
        )
        self.B_build[1] += [self.contact_count] * 7
        self.B_build[2] += list(range(7 * j, 7 * j + 7))

        self.d_overlaps[self.contact_count] = overlap_dist

        self.contact_count += 1
        # TODO: raise warning pause instead
        if self.contact_count >= self.max_contacts:
            raise RuntimeError(
                f"Reached max contact count of {self.max_contacts}. Stopping simulation!"
            )

    def build_B_filament(
        self,
        joint: CellState,
        filament: CellState,
        contact: np.ndarray,
        overlap_dist: np.floating,
    ) -> None:
        self.build_B(joint, filament, contact, overlap_dist)
        for i in [0]:
            for j in [-8, -1]:
                self.B_build[i][j] = 0
        # i.e. no growth pressure terms

    def find_contact_filament(
        self, joint: CellState, f1: CellState, f2: CellState
    ) -> tuple[np.ndarray, np.floating, np.ndarray, np.floating]:
        c1 = joint.pos
        p1, p2 = f1.ends
        if np.linalg.norm(p1 - c1) < np.linalg.norm(p2 - c1):
            c2 = p1
        else:
            c2 = p2
        p1, p2 = f2.ends
        if np.linalg.norm(p1 - c1) < np.linalg.norm(p2 - c1):
            c3 = p1
        else:
            c3 = p2

        # n1 = (c1 - f1.pos) / np.linalg.norm(c1 - f1.pos)
        # n2 = (c1 - f2.pos) / np.linalg.norm(c1 - f2.pos)
        n1 = (c1 - c2) / np.linalg.norm(c1 - c2)
        n2 = (c1 - c3) / np.linalg.norm(c1 - c3)
        # n1 += f1.dir * np.dot(f1.dir, n1)
        # n2 += f2.dir * np.dot(f2.dir, n2)

        # spring_force = 0.5 * np.cross(n1, n2) / (1 - np.dot(n1, n2)) * np.cross(n1, n2)
        # n1 += (f1.joint_growth - f1.growthRate) * np.cross(n1, -spring_force)
        # n2 += (f2.joint_growth - f2.growthRate) * np.cross(n2, spring_force)
        # n1 += np.cross(n1, -spring_force)
        # n2 += np.cross(n2, spring_force)

        # n1 = n1 / np.linalg.norm(n1)
        # n2 = n2 / np.linalg.norm(n2)

        d1 = np.linalg.norm(c2 - c1 + n1 * joint.radius)
        d2 = np.linalg.norm(c3 - c1 + n2 * joint.radius)

        x1 = -1 if np.linalg.norm(c2 - c1) < joint.radius else 1
        x2 = -1 if np.linalg.norm(c3 - c1) < joint.radius else 1
        # x1, x2 = 1, 1

        # print(joint.time, c1, c2, c3, n1, n2, d1, d2)

        c1 = c1 - c1
        c2 = c2 - f1.pos
        c3 = c3 - f2.pos

        strength = 100

        return (
            np.array([strength * n1, c1, c2]),
            strength * x1 * d1,
            np.array([strength * n2, c1, c3]),
            strength * x2 * d2,
        )
