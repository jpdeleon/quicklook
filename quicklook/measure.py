"""
python version of skimage.measure.find_contours:
https://github.com/scikit-image/scikit-image/blob/main/skimage/measure/_find_contours.py
"""

import numpy as np
from collections import deque

__all__ = ["find_contours"]


def find_contours(
    image,
    level=None,
    fully_connected="low",
    positive_orientation="low",
    *,
    mask=None,
):
    """Find iso-valued contours in a 2D array for a given level value.

    Uses the "marching squares" method to compute a the iso-valued contours of
    the input 2D array for a particular level value. Array values are linearly
    interpolated to provide better precision for the output contours.

    Parameters
    ----------
    image : 2D ndarray of double
        Input image in which to find contours.
    level : float, optional
        Value along which to find contours in the array. By default, the level
        is set to (max(image) + min(image)) / 2

        .. versionchanged:: 0.18
            This parameter is now optional.
    fully_connected : str, {'low', 'high'}
         Indicates whether array elements below the given level value are to be
         considered fully-connected (and hence elements above the value will
         only be face connected), or vice-versa. (See notes below for details.)
    positive_orientation : str, {'low', 'high'}
         Indicates whether the output contours will produce positively-oriented
         polygons around islands of low- or high-valued elements. If 'low' then
         contours will wind counter- clockwise around elements below the
         iso-value. Alternately, this means that low-valued elements are always
         on the left of the contour. (See below for details.)
    mask : 2D ndarray of bool, or None
        A boolean mask, True where we want to draw contours.
        Note that NaN values are always excluded from the considered region
        (``mask`` is set to ``False`` wherever ``array`` is ``NaN``).

    Returns
    -------
    contours : list of (n,2)-ndarrays
        Each contour is an ndarray of shape ``(n, 2)``,
        consisting of n ``(row, column)`` coordinates along the contour.
    """
    if fully_connected not in ("high", "low"):
        raise ValueError(
            'Parameters "fully_connected" must be either ' '"high" or "low".'
        )
    if positive_orientation not in ("high", "low"):
        raise ValueError(
            'Parameters "positive_orientation" must be either '
            '"high" or "low".'
        )
    if image.shape[0] < 2 or image.shape[1] < 2:
        raise ValueError("Input array must be at least 2x2.")
    if image.ndim != 2:
        raise ValueError("Only 2D arrays are supported.")
    if mask is not None:
        if mask.shape != image.shape:
            raise ValueError(
                'Parameters "array" and "mask"' " must have same shape."
            )
        if not np.can_cast(mask.dtype, bool, casting="safe"):
            raise TypeError('Parameter "mask" must be a binary array.')
        mask = mask.astype(np.uint8, copy=False)
    if level is None:
        level = (np.nanmin(image) + np.nanmax(image)) / 2.0

    segments = _get_contour_segments(
        image.astype(np.float64),
        float(level),
        fully_connected == "high",
        mask=mask,
    )
    contours = _assemble_contours(segments)
    if positive_orientation == "high":
        contours = [c[::-1] for c in contours]
    return contours


def _assemble_contours(segments):
    current_index = 0
    contours = {}
    starts = {}
    ends = {}
    for from_point, to_point in segments:
        # Ignore degenerate segments.
        # This happens when (and only when) one vertex of the square is
        # exactly the contour level, and the rest are above or below.
        # This degenerate vertex will be picked up later by neighboring
        # squares.
        if from_point == to_point:
            continue

        tail, tail_num = starts.pop(to_point, (None, None))
        head, head_num = ends.pop(from_point, (None, None))

        if tail is not None and head is not None:
            # We need to connect these two contours.
            if tail is head:
                # We need to closed a contour: add the end point
                head.append(to_point)
            else:  # tail is not head
                # We need to join two distinct contours.
                # We want to keep the first contour segment created, so that
                # the final contours are ordered left->right, top->bottom.
                if tail_num > head_num:
                    # tail was created second. Append tail to head.
                    head.extend(tail)
                    # Remove tail from the detected contours
                    contours.pop(tail_num, None)
                    # Update starts and ends
                    starts[head[0]] = (head, head_num)
                    ends[head[-1]] = (head, head_num)
                else:  # tail_num <= head_num
                    # head was created second. Prepend head to tail.
                    tail.extendleft(reversed(head))
                    # Remove head from the detected contours
                    starts.pop(head[0], None)  # head[0] can be == to_point!
                    contours.pop(head_num, None)
                    # Update starts and ends
                    starts[tail[0]] = (tail, tail_num)
                    ends[tail[-1]] = (tail, tail_num)
        elif tail is None and head is None:
            # We need to add a new contour
            new_contour = deque((from_point, to_point))
            contours[current_index] = new_contour
            starts[from_point] = (new_contour, current_index)
            ends[to_point] = (new_contour, current_index)
            current_index += 1
        elif head is None:  # tail is not None
            # tail first element is to_point: the new segment should be
            # prepended.
            tail.appendleft(from_point)
            # Update starts
            starts[from_point] = (tail, tail_num)
        else:  # tail is None and head is not None:
            # head last element is from_point: the new segment should be
            # appended
            head.append(to_point)
            # Update ends
            ends[to_point] = (head, head_num)

    return [np.array(contour) for _, contour in sorted(contours.items())]


def _get_fraction(from_value, to_value, level):
    if to_value == from_value:
        return 0
    return (level - from_value) / (to_value - from_value)


def _get_contour_segments(array, level, vertex_connect_high, mask):
    """Iterate across the given array in a marching-squares fashion,
    looking for segments that cross 'level'. If such a segment is
    found, its coordinates are added to a growing list of segments,
    which is returned by the function.  if vertex_connect_high is
    nonzero, high-values pixels are considered to be face+vertex
    connected into objects; otherwise low-valued pixels are.

    Positions where the boolean array ``mask`` is ``False`` are considered
    as not containing data.
    """

    segments = []

    use_mask = mask is not None
    square_case = 0
    # cdef tuple top, bottom, left, right
    # cdef cnp.float64_t ul, ur, ll, lr
    # cdef Py_ssize_t r0, r1, c0, c1

    for r0 in range(array.shape[0] - 1):
        for c0 in range(array.shape[1] - 1):
            r1, c1 = r0 + 1, c0 + 1

            # Skip this square if any of the four input values are masked out.
            if use_mask and not (
                mask[r0, c0] and mask[r0, c1] and mask[r1, c0] and mask[r1, c1]
            ):
                continue

            ul = array[r0, c0]
            ur = array[r0, c1]
            ll = array[r1, c0]
            lr = array[r1, c1]

            # Skip this square if any of the four input values are NaN.
            if np.isnan(ul) or np.isnan(ur) or np.isnan(ll) or np.isnan(lr):
                continue

            square_case = 0
            if ul > level:
                square_case += 1
            if ur > level:
                square_case += 2
            if ll > level:
                square_case += 4
            if lr > level:
                square_case += 8

            if square_case in [0, 15]:
                # only do anything if there's a line passing through the
                # square. Cases 0 and 15 are entirely below/above the contour.
                continue

            top = r0, c0 + _get_fraction(ul, ur, level)
            bottom = r1, c0 + _get_fraction(ll, lr, level)
            left = r0 + _get_fraction(ul, ll, level), c0
            right = r0 + _get_fraction(ur, lr, level), c1

            if square_case == 1:
                # top to left
                segments.append((top, left))
            elif square_case == 2:
                # right to top
                segments.append((right, top))
            elif square_case == 3:
                # right to left
                segments.append((right, left))
            elif square_case == 4:
                # left to bottom
                segments.append((left, bottom))
            elif square_case == 5:
                # top to bottom
                segments.append((top, bottom))
            elif square_case == 6:
                if vertex_connect_high:
                    segments.append((left, top))
                    segments.append((right, bottom))
                else:
                    segments.append((right, top))
                    segments.append((left, bottom))
            elif square_case == 7:
                # right to bottom
                segments.append((right, bottom))
            elif square_case == 8:
                # bottom to right
                segments.append((bottom, right))
            elif square_case == 9:
                if vertex_connect_high:
                    segments.append((top, right))
                    segments.append((bottom, left))
                else:
                    segments.append((top, left))
                    segments.append((bottom, right))
            elif square_case == 10:
                # bottom to top
                segments.append((bottom, top))
            elif square_case == 11:
                # bottom to left
                segments.append((bottom, left))
            elif square_case == 12:
                # lef to right
                segments.append((left, right))
            elif square_case == 13:
                # top to right
                segments.append((top, right))
            elif square_case == 14:
                # left to top
                segments.append((left, top))

    return segments
