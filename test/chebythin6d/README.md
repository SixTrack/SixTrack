# Overview of the tests
* `thin6d` is very similar to that of the elens. Test based on LHC lattice with 4 CHEBY lenses inserted at IP4. The same set of coefficients is used for all lenses, but with 3 different ways of inputting them;
* `thin6d_DYNK` is like `thin6d` but the four lenses are assigned time-varying kicks - switch off, pulsing, random on/off, random kick amplitude;
* `thin6d_ions` is very similar to `thin6d` but ions are tracked, and for more turns (i.e. 10k);
* `thin6d_ramp_DYNK` is like `thin6d_DYNK` but three (out of four) lenses are assigned time-varying kicks during a ramp;
* `thin6d_FOX_kick` is like `thin6d` but the kick computed by FOX is plotted against the kicks computed during tracking. A simpler map is used, representing a quadrupole with different gradients on the vertical and horizontal planes. Shifts and rotation angles are kept as in the original test;
* `thin6d_FOX_tune` is like `thin6d` but the test is used only to compute the tune by only a chebyshev map. The chebyshev lens implements a regular quad, which should induce a tune shift of the order of 2E-4;

# Python scripts
The python script `cheby_plot_kick.py` plots the kick received from the chebyshev leses. In total 4 lenses are inserted to check 4 different cases:
1. `cheby1`: uses the map in cheby.dat as is
1. `cheby2`: as `cheby1` with offset of (-2, 2) mm (x,y)
1. `cheby3`: as `cheby1` with offset of ( 1,-1) mm (x,y), rotated by -90 degs
1. `cheby4`: as `cheby1` rotated by 160 degs

The particle coordinates are dumped before and after the 4 lenses in the files `CHEBY_DUMP_*`.
The kick given by the lens is then just the difference between the particle coordinates after and before the elens.
Various R1 and R2 are set. For further infos, please see `fort.3`.

The script checks that x,y are unchanged (`rrin-rrout=0`) and then plots the difference in (x',y').
Offsets and tilt angle should be visible in maps.
The script also checks that points outside the domain of the lenses are not assigned a kick, and points inside the domain are always assigned a kick.

For the `DYNK` and `ion` tests, the script makes also a plot of the kicks along a given theta in the map ref frame.
The four sets of plots should be identical, apart from the time dependence.
In addition, for `ions`, the plot disentangles the different ion species.

The test also dumps the potential map of each lens - a python script is available to plot and check all of them.
