"""
@brief Contact handling for filamenting bacterial cells.
"""

import numpy as np

from CellModeller.CellState import CellState
from CellModeller.Modules.Biophysics.Contacts.BacterialContacts.BacillusContact import (
    BacillusContact,
)
from CellModeller.Modules.Biophysics.Contacts.Contact import ContactType

## @brief Relative strength of filament-joint interaction.
# @details If this value is too low ~1 then filamenting cells are easily broken
#       apart by other mechanical forces. If it is too high then these interactions
#       dominate in the solver and no other interactions take place.
FILAMENTING_STRENGTH = 1e2


class FilamentingContact(BacillusContact):
    """
    @brief Asymmetric contacts between a filament and joint.
    """

    def find_contact(self, cell: CellState) -> ContactType:
        """
        @brief What I technically spent the whole term working on...
        @details Defer to BacillusContact if not a filament, otherwise,
                determine which cell in the pair is the joint and filament and
                pass to find_contact_filament() in the correct order.
        """
        if self.cell.joint and cell.index in self.cell.filaments:
            dist, contacts = self.find_contact_filament(self.cell, cell)
            return (dist, (contacts[0], contacts[1], contacts[2]))
        elif cell.joint and self.cell.index in cell.filaments:
            dist, contacts = self.find_contact_filament(cell, self.cell)
            return (dist, (-contacts[0], contacts[2], contacts[1]))
        else:
            return super().find_contact(cell)

    def find_contact_filament(
        self, joint: CellState, filament: CellState
    ) -> ContactType:
        """
        @brief Simplified version of BacillusContact.closest_point_on_segments().
        @details Length 0 joints are single points and filaments are always anchored
                at ends, therefore there are only 2 points to check. Also, instead
                of ignoring distances further than contact range, interpret as a
                negative overlap to force a joint and segment back together. Scale
                overlap and normal by FILAMENTING_STRENGTH, cancels out in solver
                resulting in correct delta x but more weight than applied forces.
        """
        c1 = joint.pos
        p1, p2 = (
            filament.pos + 0.5 * filament.length * filament.dir,
            filament.pos - 0.5 * filament.length * filament.dir,
        )
        if np.linalg.norm(p1 - c1) < np.linalg.norm(p2 - c1):
            c2 = p1
        else:
            c2 = p2

        n = (c1 - c2) / np.linalg.norm(c1 - c2)
        d = np.linalg.norm(c2 - c1 + n * joint.radius)
        x = -1 if np.linalg.norm(c2 - c1) < joint.radius else 1

        c1 = c1 - c1
        c2 = c2 - filament.pos

        return (FILAMENTING_STRENGTH * d * x, (FILAMENTING_STRENGTH * n, c1, c2))
