"""
@brief Biophysics handling for filamenting bacterial cells.
"""

from collections.abc import Sequence
from typing import Optional

import numpy as np

from CellModeller.CellState import CellState
from CellModeller.Modules.Biophysics.BacterialModels.BacillusModel import BacillusModel
from CellModeller.Modules.Biophysics.Contacts.BacterialContacts.FilamentingContact import (
    FilamentingContact,
)
from CellModeller.Modules.Biophysics.Contacts.Contact import Contact


class FilamentingModel(BacillusModel):
    """
    @brief Extension of BacillusModel with identical underlying physics.
    @details Applies the same matrix-based method, providing different entries
            for cells marked as filamenting to replicate observed behaviour.
            In addition to those used by BacillusModel and BacillusRenderer,
            this class defines the following cell attributes:
                - filamenting
                - joint
                - filaments
                - joints
                - neighbouring_filaments
                - restoring_strength
    """

    def __init__(
        self,
        entities: list[Contact] = [],
        jitter_z: bool = False,
        bin_width: int = 5,
        atol: float = 0.000001,
        btol: float = 0.000001,
        maxiter: Optional[int] = None,
        conlim: Optional[float] = None,
        muA: float = 1,
        gamma: float = 1000000,
        alpha: float = 10000,
        jitter_magnitude: float = 0.001,
    ) -> None:
        super().__init__(
            entities,
            jitter_z,
            bin_width,
            atol,
            btol,
            maxiter,
            conlim,
            muA,
            gamma,
            alpha,
            jitter_magnitude,
        )
        add_numpy_attrs = np.dtype(
            [
                ("filamenting", "?"),
                ("joint", "?"),
                ("filaments", "2i"),
                ("joints", "2i"),
                ("neighbouring_filaments", "2i"),
                ("restoring_strength", "f"),
            ]
        )
        if self.cell_attrs is not None:
            self.cell_attrs = np.dtype(self.cell_attrs.descr + add_numpy_attrs.descr)
        self.cell_contact = FilamentingContact

    def compute_contact(
        self, cell1: CellState, cell2: CellState, cell_states: Sequence[CellState]
    ) -> None:
        """
        @brief Modified contact determination including filamenting cells.
        @details Handle 3 cases: a) regular contacts that defer to BacillusModel
                methods, b) filamenting segments and their joints that are resolved
                using FilamentingContact and do not apply a growth resistance to
                either filament (zero out 7th term of generalized coords), c) two
                filaments sharing the same joint added via build_B_joint_torsion()
                which is modified to generate a growth torque force.
        """
        if (
            cell1.filamenting
            and (not cell1.joint)
            and cell2.index in cell1.neighbouring_filaments
        ):
            joint_index = int((set(cell1.joints) & set(cell2.joints)).pop())
            self.build_B_joint_torsion(cell1, cell2, cell_states[joint_index])
        else:
            super().compute_contact(cell1, cell2, cell_states)
            if (cell1.joint and cell2.index in cell1.filaments) or (
                cell2.joint and cell1.index in cell2.filaments
            ):
                for i in [0]:
                    for j in [-8, -1]:
                        self.B_build[i][j] = 0

    def build_B_joint_torsion(self, f1: CellState, f2: CellState, joint: CellState):
        """
        @brief Compute growth torsion for two neighbouring filaments.
        @details Obtain the relative directions of each filament from the joint,
                match direction of segment in order to apply torques. Torque
                scaled by (1 + a . b) / (1 - a . b) such that maximum is when
                segments are antiparallel, half-max when segments at 90deg, and
                0 when segments are parallel. Computation of B follows similarly
                as in the method in BacillusModel, with regards to an applied
                torque at one end that displaces the center of the cell.
        """
        cells = [f1, f2]
        rel_dirs = np.array([cell.pos - joint.pos for cell in cells])
        rel_dirs = rel_dirs / np.linalg.norm(rel_dirs)
        normals = [
            cell.dir * np.sign(np.dot(cell.dir, rel_dirs[i]))
            for i, cell in enumerate(cells)
        ]
        torque = np.cross(*normals)
        spring_force = np.linalg.norm(
            (1 + np.dot(*normals)) * torque / (1 - np.dot(*normals))
        )

        for cell, normal, rel_dir, side in zip(cells, normals, rel_dirs, [1, -1]):
            i = cell.index
            self.B_build[0] += (
                list(-rel_dir - np.cross(normal, torque * side))
                + list(np.cross(torque * side, cell.dir))
                + [0]
            )
            self.B_build[1] += [self.contact_count] * 7
            self.B_build[2] += list(range(7 * i, 7 * i + 7))

        self.add_contact(joint.restoring_strength * spring_force)
