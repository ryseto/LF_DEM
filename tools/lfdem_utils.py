import numpy as np


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

    box_min = np.array(box_lim[:, 0])
    box_max = np.array(box_lim[:, 1])
    box_size = box_max - box_min

    lees_edwards_disp = lees_edwards_strain * box_size

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

    pairwise_separation = x[:, np.newaxis, :] - x
    periodize(pairwise_separation,
              box_lim,
              lees_edwards_strain,
              gradient_direction)
    return np.linalg.norm(pairwise_separation, axis=-1)
