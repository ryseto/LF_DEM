<h1> LF_DEM </h1>

A code simulating simple shear flow of dense,
overdamped suspensions of spherical particles.

 Physical interactions included:

Hydrodynamics
"Hard" contacts (with sliding friction)
Potential interaction
Brownian motion


<h2> Requirements </h2>

LF_DEM requires the sparse linear algebra software
[SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) to
be installed. Optionally, SuiteSparse can use [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for matrix ordering.

<h2> Installation </h2>

First, edit the Makefile and change the variables ```OS_version``` to
"OSX" or "Linux" depending on your environment. This variable controls
the include paths and flags to be used to compile. This has been
tested on very few machines, and is probably not generic.

If Metis is installed, also switch ```UseMetis``` to "yes". 

Once those changes to Makefile saved, you can simply compile in a terminal via:

```
$ make
```

<h2> Usage </h2>