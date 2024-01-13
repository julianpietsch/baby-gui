# BABY-GUI: user interfaces that play nice with the BABY

This repository contains user interfaces for Matlab that can be used to
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

