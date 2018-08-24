# A MAD-X/SixTrack output testing and comparison framework

This directory and all subdirectories are dedicated to provide a testing framework for MAD-X and SixTrack. This file is intended to serve as an explanation for how it is structured, meant to be used and could be expanded.

## What is the testing framework **exactly**?

It is simply a set of Python modules performing I/O and conversions as well as providing abstraction with the intent of automatizing the process of defining simulation runs that are equivalent accross MAD-X and SixTrack. It was developed with the idea of comparing tracking outputs between MAD-X and SixTrack but has been extended to also provide functionality for dealing with twiss files.

A quick summation of the scripts:
- compareSix2Mad.py: The bulk of the script. Sets up files, runs MAD-X/SixTrack and compares outputs.
- latticeConstructor.py: Abstraction for creating a lattice from Python.
- elementTest.py: Contains examples of how the framework can be used.

## An overview of how the testing framework works

## Directories in use
The current directory has three subdirectories: *scenarios/*, *work/* and *output/*. 'scenarios/' is where MAD-X/SixTrack input files are stored for use in different tests. Files belonging to the same test have the same prefix, e.g. 'test1.macro.3' and 'test1.macro.madx' would have prefix 'test1' and belong to the same test. These files are referred to as 'scenario' files and one such file can be seen below:

```
circum=%(circum)s;
lcell=20.;
f=lcell/sin(pi/5)/4;
k=1.0/f;

beam, particle=proton, energy = %(E0_GeV)s;
%(s_defs)s

seq: sequence, refer=entry, l=circum;
%(s_elems)s
endsequence;

use, sequence=seq;

select, flag=twiss, column=name, s, x,px, y,py, betx, bety;

OPTION, VERBOSE=true;

survey, file="survey.out";
twiss, file="fodo.twiss";
TRACK, ONEPASS=true, FILE="mad_output";
START, X=%(x_m)s, PX=%(PX)s, Y=%(y_m)s, PY=%(PY)s, PT=%(pt_mad)s, T=%(time_mad)s;
OBSERVE, PLACE="%(mad_track_element)s";
RUN, TURNS=%(nbr_turns)s;

sixtrack, CAVALL;	

stop;
```
A .madx scenario file.

The idea behind these files is that all the '%(...)s' are placeholders for strings that can be loaded into it via I/O and that subsequently these dummy files can be used as templates for simulation runs. By then making sure (via conversions and proper string formatting) that the formatted .madx files correspond to the same setup given by all \*.3, \*.2 etc. you can achieve simulation setups that are the same between MAD-X and SixTrack.

The 'work/' directory is where scenario files get copied to, formatted and then input to MAD-X and SixTrack. Because of this, 'work/' is during and after a testing run full of input and output files. Due to how this directory is used for extensive I/O including the copying and deletion of files, users storing files in this directory is heavily discouraged.

The 'output/' directory is reserved for output files requested by the user.

## 