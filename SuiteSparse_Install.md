# Installing SuiteSparse on HPCC machines #

Right now this HOWTO describes the installation of SuiteSparse (and optionally Metis) to run on CPUs on Andy and to run on GPUs Penzias.
The latest version of SuiteSparse is available [here](http://faculty.cse.tamu.edu/davis/suitesparse.html).

## Note on the use of Metis

Metis sometimes causes a segmentation fault, leading to a
crash of LF\_DEM. This is a Metis bug. We therefore do not recommend
the use of Metis with LF\_DEM at the moment.

## Andy

### Metis

Metis is optional for SuiteSparse.  If you want to use Metis, follow
the instructions in this section, otherwise jump directly to the
[SuiteSparse](#SuiteSparse) section. You can download Metis
[here](http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD). Get the
4.0.1 version (only metis-4.0.1 is cited by SuiteSparse doc, so for now better
stick with that one even if there are more recent versions).

Untar the Metis library, and move the
`metis-4.0/` directory in the SuiteSparse top directory:
``` $ mv
path/to/yr/lib/metis-4.0 path/to/yr/lib/SuiteSparse/ ```

By default metis uses cc as a compiler, which points to gcc on
Andy. You must change it to point to icc. Edit
`SuiteSparse/metis-4.0/Makefile.in` and define the compiler to be icc,
by changing `CC = cc` to `CC = icc`.

Metis naming of the log2
function clashes with a system library on recent Linux systems. To
avoid that, edit `SuiteSparse/metis-4.0/Lib/rename.h` and change
the last line `#define log2 __log2` to  `# #define log2 METIS__log2`.


Now we can go on and compile the whole SuiteSparse+Metis.

### SuiteSparse

First create `~/usr/lib` and `~/usr/include` if you don't have such directories.
```
$ mkdir ~/usr; mkdir ~/usr/local; mkdir ~/usr/include
```

Select the Makefile configuration file for a compilation on Linux:
```
$ cd path/to/yr/lib/SuiteSparse/SuiteSparse_config
$ cp SuiteSparse_config_linux.mk SuiteSparse_config.mk
```

Then edit this makefile `SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk` and apply the following changes:

- define `CC = icc` and `CXX = icpc` (on Andy it is NOT the default)
- uncomment the flags set to use MKL BLAS. You must obtain the following two sets of lines:
	+ `# for the MKL BLAS:`  
    `CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3 -fexceptions -fPIC -I$(MKLROOT)/include -D_GNU_SOURCE`  
	+ `# MKL`  
	  `BLAS = -Wl -mkl`  
	  `LAPACK = `
- if you want to install locally in the `~/usr/` directory, define `INSTALL_LIB` as `~/usr/lib` and `INSTALL_INCLUDE` as `~/usr/include`

If you are also installing Metis, edit `SuiteSparse/Makefile` and append the line `( $(CP) metis-4.0/libmetis.a $(INSTALL_LIB)/ && chmod 644 $(INSTALL_LIB)/libmetis.a  )` to the `install:` rule.

Then you can compile, by getting back to the main SuiteSparse directory and make:
```
$ cd ..
$ make
$ make install
```

### PBS scripts

Header of a job script file on Andy:
```
#!/bin/bash
#
#PBS -q production
#PBS -N your_job_name
#PBS -l select=1:ncpus=1
#PBS -l place=free
#PBS -V

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export CHOLMOD_OMP_NUM_THREADS=1
```



## Penzias

### SuiteSparse

First create `~/usr/lib` and `~/usr/include` if you don't have such directories.
```
$ mkdir ~/usr; mkdir ~/usr/local; mkdir ~/usr/include
```

Select the Makefile configuration file for a compilation on Linux for GPU:
```
$ cd path/to/yr/lib/SuiteSparse/SuiteSparse_config
$ cp SuiteSparse_config_GPU_icc.mk SuiteSparse_config.mk
```

Then edit this makefile and apply the following changes:

- define `CC = icc` and `CXX = icpc`
- remove the `-lrt` flag from the variable `LIB`
- if you want to install locally in the `~/usr/` directory, define `INSTALL_LIB` as `~/usr/lib` and `INSTALL_INCLUDE` as `~/usr/include`

Then you can compile, by getting back to the main SuiteSparse directory and make:
```
$ cd ..
$ make
$ make install
```

### PBS scripts

Header of a GPU job script file on Penzias:
```
#!/bin/bash
#
# My serial PBS test job.
#
#PBS -q production
#PBS -l select=1:ncpus=1:ngpus=1:mem=2880mb:accel=kepler
#PBS -N your_script_name
#PBS -l place=free
#PBS -V
```
