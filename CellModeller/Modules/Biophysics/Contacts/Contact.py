"""
@brief Namespace for classes related to biophysics module contact events.
"""

from typing import Optional, Protocol

import numpy as np

from CellModeller.CellState import CellState

## @brief Return type for find_contact() functions.
ContactType = tuple[Optional[np.floating], tuple[np.ndarray, ...]]


class Contact(Protocol):
    """
    @brief Abstract protocol for generic contacts.
    """

    def find_contact(self, cell: CellState) -> ContactType:
        """
        @brief Process a contact event with this entity.
        @param cell The contacting cell.
        @return ContactType containing contact data. Tuple of optional floating
        overlap distance of contact and a variable number of vectors related
        to the contact event. If the float is None, an event is considered not
        to be a valid contact and handled accordingly. The usual vectors returned
        are (as interpreted by BacillusModel): the contact normal on the surface
        of the cell, followed by the position of the contact relative to the
        center of the cell (cell.pos).
        """
        ...


class CellContact(Contact):
    """
    @brief Extension of Contact for cell-cell events.
    @details Calling find_contact() evaluates the contact relative to the cell
            a CellContact is initialized with. The return value will include an
            extra vector representing the contact position relative to the 2nd
            cell (passed as an argument to find_contact()). The contact normal
            for the 2nd cell will be the inverse of the initial contact normal.
    @param cell The originating cell of the event.
    """

    cell: CellState

    def __init__(self, cell: CellState) -> None:
        self.cell = cell
