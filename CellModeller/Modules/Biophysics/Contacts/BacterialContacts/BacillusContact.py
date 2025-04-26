"""
@brief Contact handling for generic bacterial models.
"""

import numpy as np

from CellModeller.CellState import CellState
from CellModeller.Modules.Biophysics.Contacts.Contact import CellContact, ContactType

## @brief Arbitrary small value used to determine if lines are parallel.
EPSILON = 1e-8


class BacillusContact(CellContact):
    """
    @brief Symmetric contacts between two rod-shaped cells.
    """

    def find_contact(self, cell: CellState) -> ContactType:
        """
        @brief Heart of the simulation <3
        @details Treat cells as line segments and then dilate by radius.
        @remark If only there was a way to make this function a little faster,
                it seems like it takes a big chunk of build_matrices()...
        """
        c1, c2, distance = closest_points_on_segments(
            *self.cell.ends, *cell.ends, epsilon=EPSILON
        )
        r1: float = self.cell.radius
        r2: float = cell.radius

        if distance >= r1 + r2:
            # Not intersecting
            return (None, ())
        else:
            # Intersecting
            normal = (
                (c2 - c1) / distance
                if distance > EPSILON
                else np.array([1.0, 0.0, 0.0])
            )
            contact_point_1 = c1 + r1 * normal
            contact_point_1 = contact_point_1 - self.cell.pos
            contact_point_2 = c2 - r2 * normal
            contact_point_2 = contact_point_2 - cell.pos
            overlap_distance = r1 + r2 - distance
            return (overlap_distance, (normal, contact_point_1, contact_point_2))


def closest_points_on_segments(
    p1, q1, p2, q2, epsilon=1e-8
) -> tuple[np.ndarray, np.ndarray, np.floating]:
    """
    @brief Find the closest points between two 3D line segments:
    @details This version is symmetric: swapping p1 and q1 (or p2 and q2) gives
            the same result. This version now also handles degenerate cases!
    @param p1 Segment 1 start.
    @param q1 Segment 1 end.
    @param p2 Segment 2 start.
    @param q2 Segment 2 end.
    @return (closest_point_on_seg1, closest_point_on_seg2, distance)
    @note This function was made with GenAI: ChatGPT using the GPT-4-turbo model,
            April 2025 release.
    """
    u = q1 - p1
    v = q2 - p2
    w0 = p1 - p2

    def closest_point_on_segment(p, a, b):
        ab = b - a
        ab_len_sq = np.dot(ab, ab)
        if ab_len_sq == 0:
            return a  # a == b
        t = np.dot(p - a, ab) / ab_len_sq
        t = np.clip(t, 0, 1)
        return a + t * ab

    # Handle degenerate cases
    if np.allclose(p1, q1):
        if np.allclose(p2, q2):
            if np.allclose(p1, p2):
                return p1, p2, np.float32(0.0)  # same point
            else:
                return (
                    p1,
                    p2,
                    np.linalg.norm(p1 - p2),
                )  # both are points, not the same
        else:
            c2 = closest_point_on_segment(p1, p2, q2)
            return p1, c2, np.linalg.norm(p1 - c2)
    elif np.allclose(p2, q2):
        c1 = closest_point_on_segment(p2, p1, q1)
        return c1, p2, np.linalg.norm(c1 - p2)
    elif np.allclose(p1, p2) or np.allclose(q1, q2):
        return p1, p1, np.float32(0.0)  # any of the equal points will do

    a = np.dot(u, u)
    b = np.dot(u, v)
    c = np.dot(v, v)
    d = np.dot(u, w0)
    e = np.dot(v, w0)

    denom = a * c - b * b

    # Initialize parameters
    s, t = 0.0, 0.0

    if denom > epsilon:
        s = (b * e - c * d) / denom
        t = (a * e - b * d) / denom

        s = np.clip(s, 0.0, 1.0)
        t = np.clip(t, 0.0, 1.0)
    else:
        # Segments are nearly parallel
        s = 0.0
        t = np.dot(v, -w0) / c if c > epsilon else 0.0
        t = np.clip(t, 0.0, 1.0)

    # Recalculate closest points after clamping
    c1 = p1 + s * u
    c2 = p2 + t * v

    # Ensure we handle edge clamping correctly
    def clamp_to_segment(p, q, x):
        t = np.dot(x - p, q - p) / np.dot(q - p, q - p)
        t = np.clip(t, 0.0, 1.0)
        return p + t * (q - p)

    # Because clamping s or t affects the other, we iterate once more
    c1 = clamp_to_segment(p1, q1, c2)
    c2 = clamp_to_segment(p2, q2, c1)

    return c1, c2, np.linalg.norm(c1 - c2)
