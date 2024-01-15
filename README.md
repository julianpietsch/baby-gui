# BABY-GUI: user interfaces that play nice with the BABY

This repository contains user interfaces for MATLAB that can be used to
process, annotate and curate cell lineages in time-lapse microscopy image data
sets. The software is general and should work with any image series, but it is
particularly powerful for images of cells in microfluidic traps or trenches.
We have used it for images of yeast in ALCATRAS devices and *E. coli* in
mother machine devices. 

The user interfaces are intended to be paired with the [Birth Annotator for
Budding Yeast (BABY)](https://github.com/julianpietsch/baby) -- an image
processing pipeline for accurate single-cell growth estimation of budding
cells from bright-field stacks. The strength of BABY is in segmenting cells
that sometimes overlap (a frequent occurence for buds and mothers in budding
yeast). However, it is a general algorithm and could be trained for
segmentation of most cell types. It has been used to segment *E. coli*, for
example. 

The user interfaces here can be used both to run the BABY algorithm and to
curate cell outlines and lineage tracks. Segmented and/or curated data sets
can be exported for time-series analysis or for training the BABY algorithm.

There are two primary Graphical User Interfaces (GUIs). The first is for
flexible specification of image data for processing and/or curation. The
second is a feature-rich curation interface for annotating and/or correcting
cell outlines, tracks and lineage information. 

Main features of the curation interface:

- Simultaneously view the same cell at multiple time points, channels and Z
  sections
- Overview plot of all cell tracks to assist in the quick identification of
  tracking errors
- Edits to track or lineage labels propagate

The BABY algorithm is described in:

[Julian MJ Pietsch, Alán F Muñoz, Diane-Yayra A Adjavon, Iseabail Farquhar,
Ivan BN Clark, Peter S Swain. (2023). Determining growth rates from
bright-field images of budding cells through identifying overlaps. eLife.
12:e79812.](https://doi.org/10.7554/eLife.79812)

If you use this software (including these user interfaces), please cite us!

## Installation

### Requirements

This software requires MATLAB 2016b or later with the Image Processing Toolbox
and Parallel Computing Toolbox.

Two contributed add-ons are also required:

- "GUI Layout Toolbox" by David Sampson, and
- "JSONLab: a toolbox to encode/decode JSON files" by Qianqian Fang

These can be installed by navigating to the Add-On Explorer in MATLAB
(accesible through the "Add-Ons" drop-down menu in the "Home" of the MATLAB
Interactive Development Environment (IDE)).

### Recommended

To support a wider range of input formats/methods it is also highly
recommended to install the [OME Bio-Formats
Toolbox](https://www.openmicroscopy.org/bio-formats/downloads/). Visit the
link to download the package. Then unzip it and copy into a location that is
available on the MATLAB path (e.g., `Documents/MATLAB/Add-Ons`).

If you want to work with data-sets located on an OMERO server, you will also
need the [OMERO MATLAB language bindings
plugin](https://www.openmicroscopy.org/omero/downloads/). For servers with
older versions of OMERO, you may need to be careful to match the plugin
version with the version of OMERO you are accessing. You can find older
versions of the plugin [here](https://downloads.openmicroscopy.org/omero/).

### Set-up

To install BABY-GUI, download this repository (either by cloning or
downloading as zip - see the "<> Code" drop-down menu at the top of this
page). Then copy the code into a location that is available on the MATLAB path
(e.g., `Documents/MATLAB/Add-Ons`). 

Alternatively, add a line to your start-up script (normally located in
`Documents/MATLAB/startup.m`; simply create it if it isn't there):

```matlab
addpath(genpath('path/to/baby-gui'));
```

## Quick start

To define a new experiment (time-lapse data-set) to work on, run
`babyCreateGUI;` at the Matlab command line.

To curate the experiment, enter `babyGUI;` at the Matlab command line, and in
the pop-up window, select the `cExperiment.mat` file that was created using
the `babyCreateGUI`.
