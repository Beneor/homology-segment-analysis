# Homolgy segment analysis

## Description
This script produces a set of short (5 - 30 by default) random fragments for the specified chromosome segment,
splits chromosome into a set of long (10 kb by default) chunks,
seeks random segments within chromosome parts,
computes frequencies of occurence for segments of different length within chromosome parts,
makes table for the normalized frequencies,
and calculates Pearson correlation between fragment frequencies for chromosome discs and ectopic contact frequencies


## Installation

Just checkout from maser branch code and run master script segmentanalysis.py

## Usage

### Default setings
`segmentanalysis.py <options> chromosomeseq segment`

The required parameters are

  `chromosomeseq` ---        TXT file containing chromosome sequence (one string,
                        small letters, no name, no spaces)
                        
  `segment` ---              segment of chromosome to analyze (16113516:16900779)
                        in BED-file notation (starting from 0, end is not
                        included)

### Optional arguments

*  `-h, --help`  ---          show this help message and exit
*  `-v, --verbose`  ---       Print additional information to stdout
*  `-s, --fragmentsizes` ---
                        Set of fragment sizes to search
*  `-d, --fragmentdensity`  ---
                        Average frequency of fragments in letters Default 25
                        meand that in average each 25-th letter will be start
                        of fragment
*  `-c, --chunk`  ---
                        Chunk size to divide chromosome
*  `-i, --iterations`  ---
                        Number of iterations for fragment splitting

## Licence
This program is licensed under GNU GPLv3 licence 
