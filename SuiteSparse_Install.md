# Installing SuiteSparse on HPCC machines #

Right now this HOWTO describes the installation of SuiteSparse (and optionally Metis) to run on CPUs on Andy and to run on GPUs Penzias.
The latest version of SuiteSparse is available [here](http://faculty.cse.tamu.edu/davis/suitesparse.html).

## Andy  ##

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
$ mkdir ~/usr; mkdir ~usr/local; mkdir ~usr/include
```

Select the Makefile configuration file for a compilation on Linux:
```
$ cd path/to/yr/lib/SuiteSparse/SuiteSparse_config
$ cp SuiteSparse_config_linux.mk SuiteSparse_config.mk
```

Then edit this makefile and apply the following changes:

- define `CC = icc` and `CXX = icpc` (on Andy it is NOT the default)
- remove the `-lrt` flag from the variable `LIB`
- uncomment the flags set to use MKL BLAS. You must obtain the following two sets of lines:
	+ `# for the MKL BLAS: <br/>  
              CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3 -fexceptions -fPIC -I$(MKLROOT)/include -D_GNU_SOURCE`
	+ `# MKL
               BLAS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
         	   LAPACK =`
- if you want to install locally in the `~/usr/` directory, define `INSTALL_LIB` as `~/usr/lib` and `INSTALL_INCLUDE` as `~/usr/include`

Then you can compile, by getting back to the main SuiteSparse directory and make:
```
$ cd ..
$ make
$ make install
```

## Penzias

First create `~/usr/lib` and `~/usr/include` if you don't have such directories.
```
$ mkdir ~/usr; mkdir ~usr/local; mkdir ~usr/include
```

Select the Makefile configuration file for a compilation on Linux:
```
$ cd path/to/yr/lib/SuiteSparse/SuiteSparse_config
$ cp SuiteSparse_config_GPU_icc.mk SuiteSparse_config.mk
```

Then edit this makefile and apply the following changes:
- define `CC = icc` and `CXX = icpc` (on Andy it is NOT the default)
- remove the `-lrt` flag from the variable `LIB`
- uncomment the flags set to use MKL BLAS. You must obtain the following two sets of lines:
  * `# for the MKL BLAS:
       CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3 -fexceptions -fPIC -I$(MKLROOT)/include -D_GNU_SOURCE`
  * `# MKL
       BLAS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
	   LAPACK =`
- if you want to install locally in the `~/usr/` directory, define `INSTALL_LIB` as `~/usr/lib` and `INSTALL_INCLUDE` as `~/usr/include`

Then you can compile, by getting back to the main SuiteSparse directory and make:
```
$ cd ..
$ make
$ make install
```

