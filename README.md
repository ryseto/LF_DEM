# LF_DEM

A code simulating simple shear flow of dense,
overdamped suspensions of spherical particles.

 Physical interactions included:

- Hydrodynamics
- "Hard" contacts (with sliding friction)
- Potential interaction OR Critical load
- Brownian motion


##Requirements

LF_DEM requires the sparse linear algebra software
[SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) to
be installed.

To install SuiteSparse on CUNY-HPPC, see [these instructions](./SuiteSparse_Install.md).

##Getting the code

### By using git

We strongly encourage to use the version control system Git (for [Mac OS X](http://git-scm.com/download/mac), for [Linux](http://git-scm.com/download/linux)). This will allow you to get updates
of the code in a very easy and clean way and to contribute to the code
by sending bug fixes or new features with a minimal effort.

You can get the code by typing in a terminal:
```
$ git clone https://rmari@bitbucket.org/rmari/lf_dem.git
```
This will download the current sources (and also the past sources).

If at any point in the future you want the latest sources of the code:
```
$ git pull
```
and that's it!

### By direct download

Download from here (bitbucket.org) by clicking the download icon (cloud), select "Download repository" and unzip.


##Installation

In the `LF_DEM` folder, edit the Makefile and change the variables
```OS_version``` to "OSX" or "Linux" depending on your
environment. This variable controls the include paths and flags to be
used to compile. This has been tested on very few machines, and is
probably not generic. If you want to install `LF_DEM` to a specific
location listed in your `$PATH`, you can set the variable
`install_dir` in the Makefile.

Once those changes to Makefile saved, you can simply compile in a terminal via:

```
$ make
$ make install
```

The second command is only needed if you want to install `LF_DEM` in the given `install_dir`.

##Usage

### Unit scales

To understand the command-line behavior of `LF_DEM`, it is useful to
know how `LF_DEM` deals with unit scales internally.

`LF_DEM` can internally work with different unit scales depending on the conditions. These unit scales are detailed in [this file](lf_dem/units.rst).


###Running a simulation

LF_DEM takes two kinds of inputs, command-line arguments and parameter
files.  The general syntax is:
```
$ LF_DEM [-p Peclet_Num ] [-c Scaled_Critical_Load ] [-r Scaled_Repulsion ] [-s Stress ] [-a Scaled_Cohesion ] [-S Stress_Sequence ] [-k kn_kt_File] [-i Provisional_Data] Configuration_File Parameter_File
```
where the only
two required arguments are the initial positions in `Configuration_File` and
the (many) simulation parameters file `Parameter_File`.

All the others are optional. If none is given, the code runs only with the hydrodynamics and the contacts (Peclet=0, no repulsion), under rate-controlled conditions (note that in this mode the code is strictly shear-rate independant). Some options are incompatible, LF_DEM should inform you so at running time.


#### Configuration file
The initial configuration can be generated by LF_DEM, see [Initial configurations](#initial).

#### Parameter file

The list of all possible parameters and their description is available
in html format in `LF_DEM/html/struct_parameter_set.html`. This file
can also be viewed online
[here](http://rmari.bitbucket.org/LF_DEM_doc/struct_parameter_set.html). (If
none of this works for you, the complete list of parameters is kept in
`LF_DEM/ParameterSet.h`.)

Although none of these parameters is
compulsory (the simulation can run with default hard-coded values), as
much as possible they should be provided by the user. One example of input parameter file is
given in the file `nobrownian_2D.txt`.

#### Rate-controlled mode

Out of the purely hydrodynamic+contact mode, there are so far three
interactions to introduce a shear-rate dependance, Brownian motion
(option `-p`), exponential repulsion (`-r`) and critical load friction
(`-c`). These command-line options have to be followed by a value. The
way LF_DEM interpret these values depend on the blend of
interactions. In Brownian mode (whenever `-p Peclet` is given), the
shear rate is given by the Peclet number. The values after `-r` or
`-c` are then interpreted as the ratio of the repulsion force or
critical load to the Brownian force scale. In non-Brownian mode, the
values after `-r` or `-c` are the non-dimensionalized shear rate (=
ratio of the hydrodynamic force and the repulsive force/critical
load).

#### Stress-controlled mode

It is selected by `-s` followed by the value of the stress. This is
simpler as so far it works only with the repulsive force. The `-r`
options should not be given in this case. `-p` and `-c` will be rejected.

#### Other options

Option                   | Role
-------------------------|--------------------------------------------
`-k  kn_kt_File`         | list of `volume_fraction kn kt dtmax` to use volume fraction dependent spring constants
`-i Provisional_Data`    | expected shear rates in stress-controlled mode to tune the output frequency
`-S Stress_Sequence`     | a sequence of `strain stress` to be followed by LF_DEM




###Initial configurations

Initial configurations can be generated through:
```
$ LF_DEM -g Random_Seed
```

LF_DEM will ask to input a series of parameters (number of particles,
dimension, etc). The generated configuration is written in a file with
a parameter dependant filename `D*N*VF*.dat`. An extra
`D*N*VF*.dat.yap` is also generated to visualize the generated
configuration with [yaplot](https://github.com/vitroid/Yaplot) or
[homer](https://github.com/rmari/homer).

##Documentation

A more complete documentation of the code is slowly building up [here](http://rmari.bitbucket.org/LF_DEM_doc/).
