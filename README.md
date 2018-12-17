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
`segmentanalysis.py <options> chromosomeseq segment [cytoamnds file]`

The required parameters are

  `chromosomeseq` ---        TXT file containing chromosome sequence (one string,
                        small letters, no name, no spaces)
                        
  `segment` ---              segment of chromosome to analyze. See below format notes

'segment' option can have flexible format to handle multiple cases. 
In general it contains name of fasta file and location of segment in this file in BED notation (chromosome, start,stop) divided by `:`
File name can be omitted to take segment from same genome as analyzed one
Start and stop locations also can be omiited to use whole fastea record. 

#### Examples
* examples/S-LIMK1.fa:S-LIMK1:500:1500 -- part of S-LIMK1 sequence from file examples/S-LIMK1.fa (letters 500-1500)
* examples/S-LIMK1.fa:S-LIMK1 -- whole S-LIMK1 sequence from file examples/S-LIMK1.fa
* :X:11982050:12772070 - 11AB -- region 11AB from Drosophila X-chromosome (same genome as analyzed)
* :X:16113516:16900779 - 14B -- different region from Drosophila X-chromosome - no correlation between homology and ectopic contacts

                        
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
