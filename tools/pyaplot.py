import numpy as np


def hfill_array(cmd_array):
    """
        Purpose:
            From an array of strings with shape (N,m), with m<=7,
            get an array of strings with shape (N,7).
    """
    fill_nb = 7 - cmd_array.shape[1]
    height = cmd_array.shape[0]
    filler = np.tile([''], (height, fill_nb))
    return np.column_stack((cmd_array, filler))


def cmd(switch_type, switch_value):
    """
        Purpose:
            Get an array of strings for commands of type switch_type taking
            values switch_value.
    """

    sval = np.array(switch_value, dtype=np.str, ndmin=1)
    switch_cmd = np.empty(sval.shape[0], dtype=np.str)
    switch_cmd[:] = switch_type
    cmd = np.column_stack((switch_cmd, sval))
    return hfill_array(cmd)


def add_cmd(yap_array, switch_type, switch_value):
    """
        Purpose:
            Append to yap_array an array of strings for commands of
            type switch_type taking values switch_value.
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
        Use case: you want to change state (e.g. width) for every object.
        You can do that in an array-like fashion, generating cmd and switch
        arrays separately, and blending them afterwards.
        This is what this function is for.
    """
    return np.reshape(np.column_stack((switch, cmd)),
                      (2*switch.shape[0], switch.shape[1]))


def interleave(array_tuple):
    """
    Purpose:

        Use case: you want to change state (e.g. width) for every object.
        You can do that in an array-like fashion, generating cmd and switch
        arrays separately, and blending them afterwards.
        This is what this function is for.

    Returns:
        [array_tuple[0][0]
         array_tuple[1][0]
         ...
         array_tuple[n][0]
         array_tuple[0][1]
         array_tuple[1][1]
         ...
         array_tuple[n][N]]
    """
    array_height = len(array_tuple)*array_tuple[0].shape[0]
    array_width = array_tuple[0].shape[1]
    return np.reshape(np.column_stack(array_tuple),
                      (array_height, array_width))


def sticks_yaparray(r1r2, thicknesses):
    """
        Get yaplot commands (as an aray of strings) to display sticks
        for each interactions with end points in r1r2.
        thicknesses contains the width of the sticks.
    """

    interaction_sticks = cmd('s', r1r2.astype(np.str))
    interaction_widths = cmd('r', thicknesses)
    yap_out = pair_cmd_and_switch(interaction_sticks, interaction_widths)

    return yap_out


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
