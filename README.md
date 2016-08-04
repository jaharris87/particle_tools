particle_tools
============

Tools for reading, analyzing, and plotting tracer particle data from 2D [CHIMERA](http://chimerasn.org) simulations

## Getting started

Particle data must be supplied by the user.

Individual ASCII files of particle histories can be generated using M.A. Chertkow's [tracer_reader](http://eagle.phys.utk.edu/chimera/trac/browser/trunk/Tools/tracer_reader)

Future versions will support HDF5 particle data


## Documentation

Some of the more complicated routines contain descritive headings and comments.


## Example

See `MATLAB/scripts/[process_B12.m](https://github.com/jaharris87/particles/blob/master/MATLAB/scripts/process_B12.m)` as an example script for setting up a working environment with particle data, loading the ASCII data into MATLAB, and initial analysis.