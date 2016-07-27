# LF_DEM cookbook

This is a set of various "how-to" for LF_DEM. They are not a complete set of examples of what you can do,
just a list of hopefully helpful cases. Not all cases will work with all versions of LF_DEM, as the code evolves constantly.
The version with which we have tested each example is given by the git commit hash at the beginning of every example. You can get to that commit with `git checkout commit_hash`.
Again, that means that this works for sure in this version, but only maybe in the previous or following versions.


### Scaling the contact model with the shear stress

**ONLY IN CONTROLLED STRESS MODE**

Tested in commit 9e2db52.

You can achieve this by setting the contact spring constants proportional to the stress, with the help of the stress unit suffix "s". A convenient way to do it is to set in the parameter file (given here in arbitrary values):
```
kn = 100s;
kt = 0.5kn;
kr = 0kn;
contact_relaxation_time = 1e-1kn;
contact_relaxation_time_tan = 0kn;
```
Note that you probably also want to scale the relaxation times with the stress, so in units of "kn" (as of this commit, you cannot give "s" units to the relaxtion times, only to the spring constants.)

When there is no other interaction than hydrodynamics and contacts, then the stress **has** to be proportional to one of the spring constants, so in that case you simply do, say, `-s 0.01kn` on the command line and launch with parameter file:
```
kt = 0.5kn;
kr = 0kn;
contact_relaxation_time = 1e-1kn;
contact_relaxation_time_tan = 0kn;
```
