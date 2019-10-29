# Tests

The four tests span different combinations of material, orientation and beam, as specified in the following list:

- `si_ch_hor_b1`: Beam 1, horizontal Si strip crystal, channeling orientation
- `si_vr_ver_b2`: Beam 2, vertical Si quasi-mosaic crystal, volume reflection orientation
- `ge_vc_ver_b1`: Beam 1, vertical Ge crystal, volume capture orientation (volume reflection still dominates, but volume capture is increased)
- `ge_am_hor_b2`: Beam 2, horizontal Ge crystal, amorphous orientation

The crystal tilt is set to different values in order to have a different dominating process in each configuration. Additionally, the geometrical parameters of the crystal are slightly different for each test.

## Python scripts

- EntCheck.py plots the particle distribution (phase-space and 2d profile) at the crystal position, before being processed by the routine
- ExitCheck.py plots the particle distribution (phase-space and 2d profile) at the crystal position, after being processed by the routine
- DeflectionCheck.py plots the deflection given to the particle by the crystal as a function of the impact angle
- InteractionCheck.py prints the percentage of particles interacting with the crystal that undergo the 4 main processes (amorphous, volume reflection, volume capture, channeling)

All the scripts are Python 3 compatible.
