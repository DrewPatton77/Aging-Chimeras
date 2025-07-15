# Changing cognitive chimera states in human brain networks with age: Variations in cognitive integration and segregation

Authors:

Drew Patton -- Code, analysis, and writer

Joern Davidsen -- Supervisor and writer

# Project Summary

The healthy brain relies on a dynamic balance between two complementary processes---integration and segregation. Integration enables coordination between distributed brain regions, while segregation reflects specialization. Optimal cognitive function emerges from a highly flexible balance between these processes as implicated by greater cognitive performance. As we age, changes in cognitive abilities occur, especially at older age, while at the same time the brain reorganizes. The exact relationship between the two is largely unknown and our work aims to tackle this challenge. Using a data-informed computational model we examine age-related differences in dynamical flexibility between segregation and integration, as captured by changes in the variable patterns of partial synchronization called chimera states---patterns where some brain regions cooperate while others remain largely independent. Our research shows that chimera states are key to maintaining the brain’s balance between integration and segregation as we age. This balance shifts differently in brain systems responsible for distinct cognitive skills, which may provide new insights into the mechanisms underlying cognitive decline and preservation in aging. In particular, it supports the idea that aging affects brain systems differently and that understanding this variability is essential for a more comprehensive view of neuro-cognitive aging. 

# Folder/File Description

* [Code](./Code)
 * [Wilson-Cowan-Simulator](./Analysis)
   This folder contains the `.m` matlab files to obtain the critical 'CE' and the synchrony matrices.   
   * [getC5crit.m](./getC5Crit.m)
     This file uses a bisection method to find the critical 'CE'.
   * [get_WC_Dynamics.m](./get_WC_Dynamics.m)
     This file runs through each individual and their single brain region stimulation and saves the resultant data.
   * [getWC.m](./getWC.m)
     This file runs the Wilson-Cowan simulation and calculates and returns the synchrony matrix.
 * [Matlab](./dir2/file21.ext)
   * [get_Patterns.m](./get_Patterns.m)
     This file runs through each synchrony matrix and applies a binarization then a Louvain algorithm to find the synchronous patterns.
   * [Tables.m](./Tables.m)
     This file runs through each structural connectome and finds its network measure then saves it as a table.
   * [Tables_Synchrony_Patterns.m](./Tables_Synchrony_Patterns.m)
     This file runs through each synchronous pattern and organizes it into a table for further analysis in `Analysis.Rmd`.
 * [Analysis](./dir2/file22.ext)
   This folder contains the rest of the analysis using the program language `R`.
   * [Analysis.Rmd](./Analysis.Rmd)
     This file is the bulk of the analysis that produces the result figures.
 * [Data](./Data)
   This folder should contain the data.
 * [Figures](./Figures)
   This folder should contain the figures once produced.
 * [Important_Files](./Important_Files)
   This folder contains supplementary file(s).
   * [Brain Region Labels 10K.txt](./Brain_Region_Labels_10K.txt)
     This text file contains the labels for each brain region and their corresponding order used in this analysis.

