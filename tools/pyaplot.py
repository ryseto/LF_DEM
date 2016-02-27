import numpy as np

def hfill_array(cmd_array):
    """
        Purpose:
            From an array of strings with shape (N,m), with m<=7,
            get an array of strings with shape (N,7).
    """
    fill_nb = 7 - cmd_array.shape[1]
    height = cmd_array.shape[0]
    filler = np.tile([''], (height,fill_nb))
    return np.column_stack((cmd_array,filler))

def cmd(switch_type, switch_value):
    """
        Purpose:
            Get an array of strings for commands of type switch_type taking values
            switch_value.
    """
    sval = np.array(switch_value, dtype = np.str,ndmin=1)
    switch_cmd = np.empty(sval.shape[0], dtype = np.str)
    switch_cmd[:] = switch_type
    cmd = np.column_stack((switch_cmd, sval))
    return hfill_array(cmd)

def add_cmd(yap_array, switch_type, switch_value):
    """
        Purpose:
            Append to yap_array an array of strings for commands of type switch_type taking values
            switch_value.
    """
    return np.row_stack((yap_array, cmd(switch_type, switch_value)))

def layer_switch(value):
    return cmd('y', value)
def add_layer_switch(yap_array, value):
    return add_cmd(yap_array, 'y', value)

def color_switch(value):
    return cmd('@', value)
def add_color_switch(yap_array, value):
    return add_cmd(yap_array, '@', value)

def radius_switch(value):
    return cmd('r', value)
def add_radius_switch(yap_array, value):
    return add_cmd(yap_array, 'r', value)


def pair_cmd_and_switch(cmd, switch):
    """
    Purpose:
        Common use case: you want to change state (e.g. width) for every object.
        You can do that in an array-like fashion, generating cmd and switch arrays separately,
        and blending them afterwards. This is what this function is for.
    """
    return np.reshape(np.column_stack((switch,cmd)),(2*switch.shape[0], switch.shape[1]))


def get_particles_yaparray(pos, rad):
    """
        Get yaplot commands (as an aray of strinfs) to display circles
        for each particle defined in (pos, rad).
        pos and rad must contain positions and radii of particles.
    """

    particle_circle_positions = cmd('c', pos)
    particle_circle_radius = cmd('r', rad)
    yap_out = pair_cmd_and_switch(particle_circle_positions, particle_circle_radius)

    return yap_out


def get_interactions_yaparray(r1r2, thicknesses):
    """
        Get yaplot commands (as an aray of strings) to display sticks
        for each interactions with end points in r1r2.
        thicknesses contains the width of the sticks.
    """

    interaction_sticks = cmd('s', r1r2.astype(np.str))
    interaction_widths = cmd('r', thicknesses)
    yap_out = pair_cmd_and_switch(interaction_sticks, interaction_widths)

    return yap_out


def get_interaction_end_points(f,p):
    """
        For each interaction in f, get the position of the particles involved.
        Positions of every particle given in p.

        Returns an array containing x1,y1,z1,x2,y2,z2 for each interaction.
    """
    # for each interaction: the particle indices
    part1 = f[:, 0].astype(np.int)
    part2 = f[:, 1].astype(np.int)

    # for each interaction: the particle positions
    r1 = p[part1, 2:5].astype(np.float)
    r2 = p[part2, 2:5].astype(np.float)

    return np.hstack((r1, r2))


def filter_interactions_crossing_PBC(f, r1r2, cutoff=4):
    """
        Exclude interactions across the boundaries.
        Return values of f and r1r2 where norm(r1r2[:,3:]-r1r2[:,:3])<cutoff.
    """
    r1 = r1r2[:, :3]
    r2 = r1r2[:, 3:]
    keep = np.linalg.norm(r2-r1, axis=1) < cutoff
    r1r2 = r1r2[keep]
    f = f[keep]
    return f, r1r2


def savetxt(outfile, yaplot_cmd_array, mode="w"):
    """
        Ouptut yaplot_cmd_array in file fname.
    """
    openedfile = True
    try:
        yap_file = open(outfile, mode+"b")
    except TypeError:
        yap_file = outfile
        openedfile = False

    np.savetxt(yap_file, yaplot_cmd_array, fmt="%s "*7)
    yap_file.write("\n".encode('utf-8'))
    if openedfile:
        yap_file.close()
