# Copyright 2016 by Tristan Bitard-Feildel.  All rights reserved.
"""MAFtk module is a set of function to be used with the MafIO parser
provided in biopython (supported in the branch MafIO of 
https://github.com/T-B-F/biopython)

This module contains read and write function to create a interval tree
based on the genomic position of the aligned blocks in MAF files.
A function is also provided to query the tree for a given species,
returning the alignments found from the query's start and stop.
The returned aligned are sliced according to the query's position.

https://github.com/T-B-F/MAFtk
"""

__version__ = "1.0"

