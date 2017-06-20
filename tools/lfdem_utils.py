import numpy as np
import dict_utils as du


def periodize(x, box_lim, lees_edwards_strain, gradient_direction=2):
    """ Periodize x according to Lees-Edwards in a box with limits box_lim.

        x is an array with shape (..., d), box_lim is an array
        with shape (d, 2) containing the min and max box coordinates
        in each direction, and lees_edwards_strain is a vector of size d.

        Note that the Lees-Edwards strain must be 0 along the
        gradient direction.
    """

    if lees_edwards_strain[gradient_direction] != 0:
        raise RuntimeError("lees_edwards_strain[gradient_direction] != 0")

    g = gradient_direction

    box_lim = np.asarray(box_lim)
    box_min = np.array(box_lim[:, 0])
    box_max = np.array(box_lim[:, 1])
    box_size = box_max - box_min

    lees_edwards_disp = np.array(lees_edwards_strain) * box_size[g]

    # start by the gradient direction, easier
    gcrossing_shift = np.array(lees_edwards_disp)
    gcrossing_shift[g] = box_size[g]

    crossing = x[..., g] > box_max[g]
    while crossing.any():
        x[crossing] -= gcrossing_shift
        crossing = x[..., g] > box_max[g]

    crossing = x[..., g] < box_min[g]
    while crossing.any():
        x[crossing] += gcrossing_shift
        crossing = x[..., g] < box_min[g]

    for d in range(x.shape[-1]):
        if d != g:
            crossing_shift = np.zeros(x.shape[-1], dtype=np.float)
            crossing_shift[d] = box_size[d]

            crossing = x[..., d] > box_max[d]
            while crossing.any():
                x[crossing] -= crossing_shift
                crossing = x[..., d] > box_max[d]

            crossing = x[..., d] < box_min[d]
            while crossing.any():
                x[crossing] += crossing_shift
                crossing = x[..., d] < box_min[d]


def pdist(x,  box_lim, lees_edwards_strain, gradient_direction=2):
    """ Get the pairwise distance matrix from an array of positions x
        according to Lees-Edwards in a box with limits box_lim.

        x is an array with shape (N, d), box_lim is an array with shape (d, 2)
        containing the min and max box coordinates in each direction,
        and lees_edwards_strain is a vector of size d.

        Note that the Lees-Edwards strain must be 0
        along the gradient direction.
    """
    box_lim = np.asarray(box_lim)
    # [-L/2,L/2] in any direction
    box_lim -= np.mean(box_lim, axis=1)[:, np.newaxis]
    lees_edwards_strain = np.asarray(lees_edwards_strain)

    pairwise_separation = x[:, np.newaxis, :] - x
    periodize(pairwise_separation,
              box_lim,
              lees_edwards_strain,
              gradient_direction)
    return np.linalg.norm(pairwise_separation, axis=-1)


def get_interaction_end_points(int_snapshot,
                               par_snapshot,
                               int_cols,
                               par_cols):
    """
        For each interaction in f, get the position of the particles involved.
        Positions of every particle given in p. Return is NOT periodized.

        Returns an array containing x1,y1,z1,x2,y2,z2 for each interaction.
    """
    p1_idx = du.matching_uniq(int_cols, ['label', '1'])[1]
    p2_idx = du.matching_uniq(int_cols, ['label', '2'])[1]
    try:
        pos_idx = du.matching_uniq(par_cols, 'position')[1]
    except Exception:
        pos_idx = slice(du.matching_uniq(par_cols, ['position', 'x'])[1],
                        du.matching_uniq(par_cols, ['position', 'z'])[1]+1)
    # for each interaction: the particle indices
    part1 = int_snapshot[:, p1_idx].astype(np.int)
    part2 = int_snapshot[:, p2_idx].astype(np.int)

    # for each interaction: the particle positions
    r1 = par_snapshot[part1, pos_idx].astype(np.float)
    r2 = par_snapshot[part2, pos_idx].astype(np.float)

    return np.hstack((r1, r2))
