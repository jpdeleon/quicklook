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
        raise ValueError('Parameters "fully_connected" must be either "high" or "low".')
    if positive_orientation not in ("high", "low"):
        raise ValueError('Parameters "positive_orientation" must be either "high" or "low".')
    if image.shape[0] < 2 or image.shape[1] < 2:
        raise ValueError("Input array must be at least 2x2.")
    if image.ndim != 2:
        raise ValueError("Only 2D arrays are supported.")
    if mask is not None:
        if mask.shape != image.shape:
            raise ValueError('Parameters "array" and "mask" must have same shape.')
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

    ul = array[:-1, :-1]
    ur = array[:-1, 1:]
    ll = array[1:, :-1]
    lr = array[1:, 1:]

    if use_mask:
        mask_ul = mask[:-1, :-1]
        mask_ur = mask[:-1, 1:]
        mask_ll = mask[1:, :-1]
        mask_lr = mask[1:, 1:]
        valid = mask_ul & mask_ur & mask_ll & mask_lr
    else:
        valid = np.ones((array.shape[0] - 1, array.shape[1] - 1), dtype=bool)

    valid &= ~(np.isnan(ul) | np.isnan(ur) | np.isnan(ll) | np.isnan(lr))

    ul = ul[valid]
    ur = ur[valid]
    ll = ll[valid]
    lr = lr[valid]

    r0, c0 = np.where(valid)
    r1 = r0 + 1
    c1 = c0 + 1

    ul = ul.ravel()
    ur = ur.ravel()
    ll = ll.ravel()
    lr = lr.ravel()

    square_case = np.zeros(len(ul), dtype=np.int32)
    square_case[ul > level] += 1
    square_case[ur > level] += 2
    square_case[ll > level] += 4
    square_case[lr > level] += 8

    # Filter out trivial cases (0 and 15) up front
    active = (square_case != 0) & (square_case != 15)
    if not active.any():
        return segments

    sc = square_case[active]
    r0_a = r0[active].astype(float)
    c0_a = c0[active].astype(float)
    r1_a = r1[active].astype(float)
    c1_a = c1[active].astype(float)
    ul_a = ul[active]
    ur_a = ur[active]
    ll_a = ll[active]
    lr_a = lr[active]

    # Vectorized fraction computation (safe divide, 0 where denom is 0)
    def _frac(fm, to):
        denom = to - fm
        return np.where(denom != 0, (level - fm) / denom, 0.0)

    top_r = r0_a
    top_c = c0_a + _frac(ul_a, ur_a)
    bot_r = r1_a
    bot_c = c0_a + _frac(ll_a, lr_a)
    lft_r = r0_a + _frac(ul_a, ll_a)
    lft_c = c0_a
    rgt_r = r0_a + _frac(ur_a, lr_a)
    rgt_c = c1_a

    # Segment lookup: for each case, (from_row, from_col, to_row, to_col)
    # Cases 6 and 9 produce two segments and depend on vertex_connect_high
    _T, _B, _L, _R = 0, 1, 2, 3  # point indices
    # Single-segment cases: case -> (from_point, to_point)
    single = {
        1: (_T, _L),
        2: (_R, _T),
        3: (_R, _L),
        4: (_L, _B),
        5: (_T, _B),
        7: (_R, _B),
        8: (_B, _R),
        10: (_B, _T),
        11: (_B, _L),
        12: (_L, _R),
        13: (_T, _R),
        14: (_L, _T),
    }
    # Rows/cols for each point type
    rows = [top_r, bot_r, lft_r, rgt_r]
    cols = [top_c, bot_c, lft_c, rgt_c]

    for case_val, (fp, tp) in single.items():
        idx = np.where(sc == case_val)[0]
        for i in idx:
            segments.append(((rows[fp][i], cols[fp][i]), (rows[tp][i], cols[tp][i])))

    # Ambiguous cases 6 and 9 (two segments each)
    if vertex_connect_high:
        pairs_6 = [(_L, _T), (_R, _B)]
        pairs_9 = [(_T, _R), (_B, _L)]
    else:
        pairs_6 = [(_R, _T), (_L, _B)]
        pairs_9 = [(_T, _L), (_B, _R)]

    for case_val, pairs in [(6, pairs_6), (9, pairs_9)]:
        idx = np.where(sc == case_val)[0]
        for i in idx:
            for fp, tp in pairs:
                segments.append(((rows[fp][i], cols[fp][i]), (rows[tp][i], cols[tp][i])))

    return segments
