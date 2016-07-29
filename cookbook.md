# LF_DEM cookbook

This is a set of various "how-to" for LF_DEM. They are not a complete set of examples of what you can do,
just a list of hopefully helpful cases. Not all cases will work with all versions of LF_DEM, as the code evolves constantly.
The version with which we have tested each example is given by the git commit hash at the beginning of every example. You can get to that commit with `git checkout commit_hash`.
Again, that means that this works for sure in this version, but only maybe in the previous or following versions.


### Scaling the contact model with the shear stress

**ONLY IN CONTROLLED STRESS MODE**

Tested in commit c87fa13.

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

### Simulation with contact dashpots and without lubrication

Tested in commit a84707e.

You can achieve this by setting `lubrication_model = none;`. In this case you have to set
both normal and tangential contact relaxation times explicitly,
because the dashpot resistance cannot be set to the resistance of lubrication at contact,
as is the default for the tangential dashpot when there is lubrication.

For instance:
```
kn = 300h;
kt = 0.5kn;
contact_relaxation_time = 4e-1kn;
contact_relaxation_time_tan = 4e-1kn;
```
Note that the value you choose for the relaxation times is important to get a physically reasonable behavior
and should be decided based on the kind of problem studied. When doing simple shear under rate controlled conditions,
an often reasonable choice (but not always, e.g. not a very low shear rates) is to give the relaxation times corresponding to small strains, like for instance:
```
contact_relaxation_time = 1e-3h;
contact_relaxation_time_tan = 2e-3h;
```

You do not have to worry too much about `lub_max_gap`. You do not have to set it to 0 for this simulation, although temporarily it is a good idea to do so for simulations with hydro and contacts only, because the range of the interactions is still set to `lub_max_gap` (such that there are a lot of useless interaction objects created if `lub_max_gap` is large, which makes the simulation needlessly slower).
