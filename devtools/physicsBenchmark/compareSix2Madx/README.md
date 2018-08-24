# A MAD-X/SixTrack output testing and comparison framework

This directory and all subdirectories are dedicated to provide a testing framework for MAD-X and SixTrack. This file is intended to serve as a brief explanation.

## What is the testing framework **exactly**?

It is simply a set of Python modules performing I/O and conversions as well as providing abstraction with the intent of automatizing the process of defining simulation runs that are equivalent accross MAD-X and SixTrack. It was developed with the idea of comparing tracking outputs between MAD-X and SixTrack but has been extended to also provide functionality for dealing with twiss files.

A quick summation of the scripts:
- compareSix2Mad.py: The bulk of the script. Sets up files, runs MAD-X/SixTrack and compares outputs.
- latticeConstructor.py: Abstraction for creating a lattice from Python.
- elementTest.py: Contains an example of how the framework can be used.

## The code should for the most part be self-explanatory

Granted time and interest, this introduction could be expanded..

J. Andersson, 24th of August 2018