# Evaluating the ambush hypothesis

### Overview

This repository contains Python code that I wrote for my BSc thesis: Evaluating the ambush hypothesis using data from 744 bacterial genomes. This worked was completed at the University of Bath in 2017.


### Abstract 

Ribosomal frameshifting events can have catastrophic consequences on the genome. The ambush hypothesis postulates that off-frame stop codons (OSCs) are selected for in the genome to catch ribosomal slippages. These ‘hidden’ stops are especially prevalent in prokaryotes, which are known to utilise programmed frameshifts and in turn may be more susceptible to random ones. Previous studies have found associations between OSCs and: GC content, codon usage, and optimal growth temperature. They have also shown variable frequencies of +1 and +2 frameshifted stops and the significant excess of OSCs in prokaryotic genomes.

The primary aims of this research were to determine the frequency with which OSCs are found in excess in bacterial genomes and to establish whether the variation of OSCs can be fully explained by the variation of GC content. The secondary aims were to discern the role of the TGA stop codon in OSC selection, assess the relative frequencies of +1 and +2 OSCs and ascertain the consistency of nucleotide intervals between OSCs.

To achieve these aims, in this study OSCs were assessed across 744 bacterial genomes split into two distinct datasets according to the stop codon permutations that can occur during translation. Levels of OSC excess were determined using a Monte Carlo approach. A significant strong negative correlation between OSC percentage and GC content was found as was a corresponding positive correlation between OSC intervals and GC content. Contrary to previous research, the majority of genomes did not show a significant excess of OSCs. Weak correlations were found between OSC excess and GC content. The TGA stop codon was determined to be the least integral to overall selection for OSCs. +2 OSCs were found to be significantly more frequent than +1 in some cases.

This research presents arguments for and against the ambush hypothesis and highlights the fact that more work is required to understand the extent and intricacies of this hypothesis.


### Scripts

[**OSC_analysis.py**](OSC_analysis.py) contains functions to:

1) Read in embl genome files and extract CDS regions
2) Check that sequences are true protein CDSs
3) Calculate Off-frame stop codon (OSC) percentage
4) Perform Monte Carlo randomisation for OSC excess calculations
5) Calculate GC and GC3 contents
6) Calculate OSC intervals

[**Execution_script.py**](Execution_script.py) contains code to excute functions in [**OSC_analysis.py**](OSC_analysis.py) and write out to a results file.
[**OSC_analysis_visualisations.R**](OSC_analysis_visualisations.R) contains code to create all of the graphs made from the output of the OSC Analysis.
