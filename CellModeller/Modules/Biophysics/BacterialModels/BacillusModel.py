"""
@brief Biophysics handling for generic bacterial models.
"""

from collections.abc import Sequence
from time import perf_counter
from typing import Any, Optional

import numpy as np
import scipy.sparse as spa
from scipy.sparse.linalg import lsmr

from CellModeller.CellState import CellState
from CellModeller.Modules.Biophysics.Contacts.BacterialContacts.BacillusContact import (
    BacillusContact,
)
from CellModeller.Modules.Biophysics.Contacts.Contact import CellContact, Contact
from CellModeller.Modules.CMModule import CMModule, cell_event


class BacillusModel(CMModule):
    """
    @brief Default biophysics module for simulating cell interactions.
    @details Implements methods for building matrices according to the method
            outlined in the original paper and included math.pdf documentation.
            Cell attributes defined by this class can be accessed by all other
            modules:
                - growth_rate
                - neighbourhood
                - ends
    @note In addition to keeping defined attributes up to date, the biophysics
            modules are also responsible for updating the following, which must
            be defined in another module or (usually) renderer to work:
                - radius
                - length
                - pos
                - dir
    @remark dir is expected to be a unit vector, although it is normalized in
            the update step of the module.
    @param entities List of non-cell entities that are a part of the simulation.
            Currently no neighbourhoods calculation like for cell so assume
            there are comparatively fewer entities in a simulation.
    @param jitter_z Whether or not to move cells in the z-direction. If set to
            False, cells will remain in the xy plane only.
    @param bin_width For calculating adjacent cells. Should be 1.25x the max
            cell length expected for the simulation. For more details, see
            Algorithms for Particle-Field Simulations with Collisions
            Hersir Sigurgeirsson, Andrew Stuart, and Wing-Lok Wan
            https://cs.uwaterloo.ca/~jwlwan/papers/SigETAL01.pdf
    @param atol Tolerance for solver Ax term. See scipy.sparse.linalg.lsmr docs.
    @param btol Tolerance for solver b term. See scipy.sparse.linalg.lsmr docs.
    @param maxiter Max solver iterations. See scipy.sparse.linalg.lsmr docs.
    @param conlim Condition limit. See scipy.sparse.linalg.lsmr docs.
    @param muA Simulation viscosity parameter. See math.pdf docs.
    @param gamma Cell growth resistance parameter. See math.pdf docs.
    @param alpha Simulation elasticity parameter. See math.pdf docs.
    @param epsilon Small value for comparing parallel lines.
    @param jitter_magnitude Amount of randomness to introduce when cells divide.
    """

    def __init__(
        self,
        entities: list[Contact] = [],
        jitter_z: bool = False,
        bin_width: int = 5,
        atol: float = 1e-6,
        btol: float = 1e-6,
        maxiter: Optional[int] = None,
        conlim: Optional[float] = None,
        muA: float = 1.0,
        gamma: float = 1e6,
        alpha: float = 1e4,
        jitter_magnitude: float = 0.001,
    ) -> None:
        self.cell_attrs = np.dtype(
            [
                ("growth_rate", "f"),
                ("neighbourhood", "3i"),
                ("ends", "2,3f"),
            ]
        )

        self.entities = entities
        self.jitter_z = jitter_z
        self.muA = muA
        self.gamma = gamma
        self.alpha = alpha
        self.jitter_magnitude = jitter_magnitude

        self.solver_args: dict[str, Any] = {"atol": atol, "btol": btol}
        if maxiter is not None:
            self.solver_args["maxiter"] = maxiter
        if conlim is not None:
            self.solver_args["conlim"] = conlim

        self.bin_width: int = bin_width
        self.neighbourhoods: dict[tuple[int, int, int], list[int]] = {}
        adjacent = (-bin_width, 0, bin_width)
        self.adjacent_bins = [
            (dx, dy, dz) for dx in adjacent for dy in adjacent for dz in adjacent
        ]
        self.contact_count: int = 0
        self.cell_contact: type[CellContact] = BacillusContact

        self.build_times: list[float] = []
        self.solver_residuals: list[float] = []
        self.A_condition: list[float] = []
        self.solver_times: list[float] = []

    def on_sim_ready(self) -> None:
        """
        @brief Allocate memory for input and output vectors of lsmr solver.
        """
        self.d_overlaps = np.zeros(self.sim_params.max_contacts, dtype=np.float32)
        self.delta_p = np.zeros(7 * self.sim_params.max_cells, dtype=np.float32)
        self.reset_matrices()

    def reset_matrices(self) -> None:
        """
        @brief Build solver input matrices in sparse COO format.
        @details Using python lists as items are only read once when converting
                so no benefit to allocating and zeroing out memory each step.
        """
        self.M_build: list[list] = [[], [], []]
        self.B_build: list[list] = [[], [], []]
        self.d_overlaps[: self.contact_count] = 0
        self.contact_count = 0

    def sim_step(self, dt: float, cell_states: Sequence[CellState]) -> None:
        """
        @brief Compute change in cell positions and update them.
        @details This occurs in 3 steps: build the contact (B) and inertia (M)
                matrices and overlaps vector (d) as outlined in math.pdf, convert
                the matrices to a sparse format and solve using lsmr iterative
                method, then push updates back to cell arrays.
        """
        count = len(cell_states)
        if count == 0:
            return

        self.build_times += [-perf_counter()]
        self.build_matrices(cell_states)
        self.build_times[-1] += perf_counter()
        if self.sim_params.verbosity_level(3):
            print(f"buildmatrices took {np.mean(self.build_times):.2f}s (avg)")
            self.build_times = []

        self.solver_times += [-perf_counter()]
        if self.contact_count > 0:
            M = spa.coo_matrix(
                (self.M_build[0], (self.M_build[1], self.M_build[2])),
                dtype=np.float32,
                shape=(7 * count, 7 * count),
            ).tocsr()
            B = spa.coo_matrix(
                (self.B_build[0], (self.B_build[1], self.B_build[2])),
                dtype=np.float32,
                shape=(self.contact_count, 7 * count),
            ).tocsr()
            d = self.d_overlaps[: self.contact_count]

            A = self.alpha * B.T @ B + M / dt
            b = -self.alpha * B.T @ d
            if self.sim_params.verbosity_level(4):
                print(f"Sparsity of A: {100*(1-A.nnz/np.prod(A.shape)):.2f}%")

            res = lsmr(A, b, **self.solver_args)
            self.delta_p[: 7 * count] = res[0]
            self.solver_residuals += [res[3]]
            self.A_condition += [res[6]]
            if self.sim_params.verbosity_level(4):
                print(f"Iterations: {res[2]}")
                print(f"Residual: {np.mean(self.solver_residuals)} (avg)")
                print(f"Condition number of A: {np.mean(self.A_condition)} (avg)")
                self.solver_residuals = []
                self.A_condition = []
        else:
            self.delta_p[: 7 * count] = 0
        self.solver_times[-1] += perf_counter()
        if self.sim_params.verbosity_level(3):
            print(f"solving took {np.mean(self.solver_times):.2f}s (avg)")

        for state in cell_states:
            self.updateCell(state, dt)

        self.reset_matrices()

    def build_matrices(self, cell_states: Sequence[CellState]) -> None:
        """
        @brief Find cell contacts and build B and M matrices.
        @details Calls build_M() on each cell, then iterates through cells
                in its neighbourhood and calls compute_contact() on pairs,
                finally check for contacts with non-cell simulation entities.
                The latter two contact checks update B via build_B().
        """
        for cell in cell_states:
            self.build_M(cell)

            for neighbour_index in self.neighbourhoods[tuple(cell.neighbourhood)]:
                if cell.index < neighbour_index:
                    neighbour = cell_states[neighbour_index]
                    self.compute_contact(cell, neighbour, cell_states)

            for entity in self.entities:
                dist, contacts = entity.find_contact(cell)
                if dist is not None:
                    normal, contact = contacts
                    self.build_B(cell, contact, normal)
                    self.add_contact(dist)

    def compute_contact(
        self, cell1: CellState, cell2: CellState, cell_states: Sequence[CellState]
    ) -> None:
        """
        @brief Determine if two cells are in contact and update B.
        @details Query the class CellContact object for actual calculation.
        """
        del cell_states
        dist, contacts = self.cell_contact(cell1).find_contact(cell2)
        if dist is not None:
            normal, contact_1, contact_2 = contacts
            self.build_B(cell1, contact_1, normal)
            self.build_B(cell2, contact_2, -normal)
            self.add_contact(dist)

    def build_M(self, cell: CellState) -> None:
        """
        @brief Add an entry to the inertia matrix.
        @details Use cell length and muA (viscosity), as well as cell dir
                to construct projection matrix for rotational component.
        """
        i = cell.index

        self.M_build[0] += [self.muA * cell.length] * 3
        self.M_build[1] += list(range(7 * i, 7 * i + 3))
        self.M_build[2] += list(range(7 * i, 7 * i + 3))

        # ||a x r||^2 = (1 - aa^T) . r (??)
        self.M_build[0] += list(
            (np.eye(3) - np.outer(cell.dir, cell.dir)).flatten()
            * self.muA
            * cell.length
            / 12
        )
        self.M_build[1] += [7 * i + 3 + a for a in [0, 0, 0, 1, 1, 1, 2, 2, 2]]
        self.M_build[2] += [7 * i + 3 + a for a in [0, 1, 2, 0, 1, 2, 0, 1, 2]]

        if cell.growth_rate > 0:
            self.M_build[0] += [self.gamma]
        else:
            self.M_build[0] += [0]
        self.M_build[1] += [7 * i + 6]
        self.M_build[2] += [7 * i + 6]

    def build_B(self, cell: CellState, contact: np.ndarray, normal: np.ndarray) -> None:
        """
        @brief Add an entry to the contact matrix.
        @details For cell-cell contacts this is only half of the entry.
        @param cell The cell in contact.
        @param contact Position of contact relative to center of cell.
        @param normal Vector of contact on surface of cell.
        """
        i = cell.index
        self.B_build[0] += (
            list(normal)
            + list(np.cross(-np.cross(normal, contact), cell.dir))
            + [np.dot(normal, cell.dir) * np.dot(contact, cell.dir)]
        )
        self.B_build[1] += [self.contact_count] * 7
        self.B_build[2] += list(range(7 * i, 7 * i + 7))

    def add_contact(self, overlap_dist: np.floating):
        """
        @brief Add an entry to the overlaps vector.
        @details Required to be called after every entry to the contact matrix,
                as its length has to match the rows of B.
        """
        self.d_overlaps[self.contact_count] = overlap_dist
        self.contact_count += 1
        if self.contact_count >= self.sim_params.max_contacts:
            self.pause()
            print(f"Reached max contact count of {self.sim_params.max_contacts}!")

    def updateCell(self, cell: CellState, dt: float) -> None:
        """
        @brief Perform cell growth and update ends.
        """
        i = cell.index

        d_pos = self.delta_p[7 * i : 7 * i + 3]
        d_dir = self.delta_p[7 * i + 3 : 7 * i + 6]
        d_len = self.delta_p[7 * i + 6]

        cell.pos += d_pos
        cell.dir += d_dir
        cell.dir = cell.dir / np.linalg.norm(cell.dir)
        growth = dt * cell.length * (cell.growth_rate / 60)
        if d_len > -growth:
            cell.length += d_len + growth

        cell.ends = (
            cell.pos + 0.5 * cell.length * cell.dir,
            cell.pos - 0.5 * cell.length * cell.dir,
        )

        # check if moved out of neighbourhood
        local_pos = cell.pos - cell.neighbourhood
        # +25% margin on all sides of bin to avoid chattering
        if any(abs(a) > 0.625 * self.bin_width for a in local_pos):
            self._unset_neighbourhood(cell)
            self._set_neighbourhood(cell)

    def new_cell(self, cell: CellState) -> None:
        """
        @brief Set default values for cell attributes used in this module.
        """
        cell.radius = 0.5
        cell.length = 3.0
        cell.dir = np.array((1, 0, 0))
        cell.growth_rate = 0.05
        cell.ends = (
            cell.pos + 0.5 * cell.length * cell.dir,
            cell.pos - 0.5 * cell.length * cell.dir,
        )
        self._set_neighbourhood(cell)

    @cell_event("divide_flag", 1)
    def divide(self, cells: tuple[CellState, ...]) -> None:
        """
        @brief Handle division events triggered by divide_flag.
        @details Also assigns priority 1 to divide so that it occurs before
                cell removal (default PreEvent of divide_flag marks parent for
                removal).
        @param cells tuple of (parent, daughter1, daughter2)
        """
        parent, *daughters = cells

        for daughter, side in zip(daughters, [1, -1]):
            daughter.radius = parent.radius
            daughter.length = parent.length / 2.0 - parent.radius
            daughter.pos = (
                parent.pos
                + (daughter.length / 2.0 + daughter.radius) * side * parent.dir
            )
            jitter = np.random.uniform(-self.jitter_magnitude, self.jitter_magnitude, 3)
            if not self.jitter_z:
                jitter[2] = 0
            daughter.dir = parent.dir + jitter
            daughter.growth_rate = parent.growth_rate
            daughter.ends = (
                daughter.pos + 0.5 * daughter.length * daughter.dir,
                daughter.pos - 0.5 * daughter.length * daughter.dir,
            )
            self._unset_neighbourhood(daughter)
            self._set_neighbourhood(daughter)

    @cell_event("remove_flag")
    def del_cell(self, cells: tuple[CellState, ...]) -> None:
        """
        @brief Handle cell removal due to remove_flag.
        """
        for cell in cells:
            self._unset_neighbourhood(cell)

    def _set_neighbourhood(self, cell: CellState) -> None:
        """
        @brief Update internal map of neighbourhoods.
        @details Called for new cells and when cells move neighbourhoods.
        """
        x, y, z = map(int, np.round(cell.pos / self.bin_width) * self.bin_width)
        cell.neighbourhood = (x, y, z)

        for dx, dy, dz in self.adjacent_bins:
            nx, ny, nz = x + dx, y + dy, z + dz
            if (nx, ny, nz) in self.neighbourhoods:
                self.neighbourhoods[(nx, ny, nz)] += [cell.index]
            else:
                self.neighbourhoods[(nx, ny, nz)] = [cell.index]

    def _unset_neighbourhood(self, cell: CellState) -> None:
        """
        @brief Remove cell from neighbourhoods.
        @details Called when removing cells and when cells move neighbourhoods.
        """
        x, y, z = cell.neighbourhood

        for dx, dy, dz in self.adjacent_bins:
            nx, ny, nz = x + dx, y + dy, z + dz
            self.neighbourhoods[(nx, ny, nz)].remove(cell.index)
            if len(self.neighbourhoods[(nx, ny, nz)]) == 0:
                del self.neighbourhoods[(nx, ny, nz)]
