
# Homology segment analysis

The application calculates sequence homology level using short fragment matches 
and identifies correlations between frequency of matched fragments (FMF) and any other genome data (for example frequency of ectopic contacts, FEC) 

The software capabilities
* Calculate homology level between provided sequence with a genome
* Calculate a correlation between homology level any numeric data specified for genome (for example FEC)
* Plot homology level graph 
* Analyse full-chromosomee correlations between FMF and any other genome data (for example FEC)  

## Licence
This program is licensed under GNU GPLv3 licence: https://www.gnu.org/licenses/gpl.html


# Installation 
Software is written using Python language version 3.
In addition to main Python package, it uses following Python modules:
* Aho-Corasick for fast fragment search
* Numpy and Scipy for correlation analysis
* Matplotlib for plotting graphs.

Both for Windows and Linux users we recommend installing all dependencies via Miniconda. The same way should work for Mac also. 
However, the application was not tested on Mac.

## Installing Anaconda
* Download installer for latest python3 version from miniconda site: https://conda.io/en/latest/miniconda.html 
* Follow the installation instructions

## Creating Conda environment
Execute following commands from your shell

    conda create --name py3 python=3.8
    conda activate py3
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda update --all
    conda install numpy scipy matplotlib pandas pyahocorasick

### Windows users note
If installation of pyahocorasick module fails, try to install Visual Studio build tools from here: https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017

### Linux users note
All Linux distributions contains Python3 by default. 
So Linux users instead of Anaconda package can install all required packages 
via distribution package manager and python pip.

Example for Ubuntu/Debian:

    sudo apt install build-essential python3-dev
    pip3 install numpy scipy matplotlib pandas pyahocorasick


## Main program installation
Just checkout from master branch.

`git clone https://bitbucket.org/beneor/homology-segment-analysis.git`

Then navigate to program base folder, and start using.

# Usage

## Calculate homology between provided sequence and genome using segmentanalysis.py

The `segmentanalysis.py` script splits segment into set of short (20 - 40 kb) fragments,
then splits chromosome into a set of long (10 kb by default) chunks, and seeks for positions where these fragments match with the genome.
Then it computes frequencies of occurrence for segments within a chromosome,
makes table of the match frequencies normalized by chunks, and calculates fragments matching frequency by cytobands, 
if cytoband information is provided.

Resulting data are grouped by genome intervals (chunks) and by cytobands (if cytoband data is available)

General run syntax is following: 

`segmentanalysis.py <options> chromosomeseq segment [cytobands]`

### Main arguments
The main arguments are

+  `fastaFileName` -      FASTA file (may be gzipped) containing genome sequence                      
+  `segment` -              segment of chromosome to analyze: file and location. See format notes below.
+  `cytobands` - optional BED-format file containing locations of cytobands or any other intervals to group homology data. 

#### Segment notation
**segment** option has a flexible format to handle multiple cases. 
In general, it contains the name of the FASTA file and location of segment in this file in BED 
notation (chromosome, start,stop) divided by **:**
File name can be omitted to take the segment from same genome as analyzed one.
Start and stop locations also can be omitted to use whole FASTA record. 

##### Examples
* `examples/S-LIMK1.fa:S-LIMK1:500:1500` -- part of S-LIMK1 sequence from file examples/S-LIMK1.fa (letters 500-1500)
* `examples/S-LIMK1.fa:S-LIMK1` -- whole S-LIMK1 sequence from file examples/S-LIMK1.fa
* `:X:11982050:12772070` - 11AB -- region 11AB from Drosophila X-chromosome (same genome as analyzed)
* `:X:16113516:16900779` - 14B -- different region from Drosophila X-chromosome which doesn't have any correlation between homology and ectopic contacts
                       
### Optional arguments
* `-h, --help`  -          Show help message and exit
* `-v, --verbose`  -       Print additional information to stdout
* `-s, --fragmentsizes` - Set of fragment sizes to search divided by comma. Default: 20,25,30,35,40
* `-d, --fragmentdensity`  - Use random chosen fragments instead of all possible.
                        Sets average frequency of fragments in letters. For
                        example 5 means that each 5-th letter will be start of
                        fragment
* `-c, --chunk`  - Chunk size to divide chromosome, in kilobases. Default: 10. 
                        
* `--nodump` - Do not save found fragments and positions to files.
* `--mindumpsize` - Minimum size of fragment to store exact locations in
                        file. Default 20


## Homology graph plotting using genomecorrelation.py 

Second script `genomecorrelation.py` plots homology level graphs (both grouped by chunks or by cytobands),
and calculated correlation of homology level with experimental data, if provided.

General run syntax is following: 

`genomecorrelation.py <options> homologyFileName [ectopicsFileName]` 

The main arguments are

+ `homologyFileName` - BED-formatted file holding homology data
+ `ectopicsFileName` - optional tab-delimeted file containing ectopic contacts data

### Optional arguments
+ `-h, --help` -           Show help message and exit
+ `-r,  --correlation` - Method for correlation calculations. Valid options are *Spearman* or *Pearson*. 
                        Default: Spearman
+ `-c, --chromosomes` - List of chromosomes divided by comma to view (view all if
                        not set). Example: 2L,X
+ `-l, --labels` -         Put labels from homology data on X axis
+ `-m, --mark` -           Add data points marks to graph
+ `-x, --xinterval` - Interval to put ticks on X-axis, in kilobases
+ `-s, --savefig` -   File name to save resulting picture

## Full-chromosome correlation analysis using analysis.py 

The script `analysis.py` analyses full-chromosomee correlations between FMF and any other genome data (for example FEC).
It includes Spearman rho correlation (R) for specific and unspecific matches, combining data and graphs plotting.

The part of the application is currently not refactored, so using it is not so convenient as other scripts
Additional_scripts
General run syntax is following: 
`analysis.py analysis.py <cytocountsData> <Xectopicsdata> [--nearest]` 

The main arguments are

+ `cytocountsData` - Folder containing pre-calculated data with cytocounts.
+ `Xectopicsdata` - optional tab-delimeted file containing ectopic contacts data

### Optional arguments
+ `strain1` - Name of the first strain to analyse. Default: CS
+ `strain2` - Name of the second strain to analyse. Default: agn
+ `-h, --help` -     Show help message and exit
+ `-s, --fragmentsizes` -- Set of fragment sizes to search. Default: CS,agn
+ `-N,  --nearest` - Perform nearest calculations. Default: No

### Data preparation
#### Cytocounts FMF data
The pre-calculated matrix of FMF should be built manually before running `analysis.py` script. 
To do this, user should first pre-calculate correlations fo different regions of the chromosome using bash script
Example script `batch-segmentanalysis-all-Xdiscs` is located in `Additional_scripts` folder.

The resulting data folder should have the following structure:
    
    cytocountsforler
        - disks
            - all
                -l<length>-cytocounts
            - nr
                -l<length>-cytocounts
        - zones
            - all
                -l<length>-cytocounts
            - nr
                -l<length>-cytocounts

The `all` and `nr` folders contain data for all fragment  and only for non-repeatable segments

The file names in the folder is consist from the zone number and disk number.
For _Drosophila melanogaster_ X-chromosome there are 19 zones and 6 disks (A-F) in each zone. 
So the file 123.txt corresponds to the disk 12C

The pre-calculated archive for _Drosophila melanogaster_ X-chromosome is available here: 
https://drive.google.com/file/d/1AOrQ06x6_oek3sfDhYZNr-pw_sNdpPe0/view?usp=sharing

#### Ectopics data
Currently the ectopic data are organized in set of 1-dimensional file tables
The data format change to more convenient 2D-matrix is under development now. 

The ectopic data folder should have the following structure:

    ectopicsdataforler
        - strain 1
            - disks
            - zones
        - strain 2
            - disks
            - zones

The file names follows the same structure as for cytocounts FMF data
* The file name for zones is just a number of zone
* The file name for disk is consist for zone (1-19) and number of disk 1-6. So the file 123.csv represents the disk 12C

The sample data archive for _Drosophila melanogaster_ X-chromosome is available here: 
https://drive.google.com/file/d/1cPDRCaTrUlxa4Pp0fR4DxnP9NzEKLCgt/view?usp=sharing

## Test run examples
Generate homology data for Drosophila genome, included in examples data set 
`./segmentanalysis.py -v -s 30 -d 2.5 examples/dm6.nounmapped.fa.gz :X:11982050:12772070 examples/DmelMapTable.160615c.bed`

View homology graph together with ectopics contacts frequency
`/genomecorrelation.py -m -c X,2L examples/dm6.nounmapped.fa.gz.X-11982050-12772070/cytocounts.l30-merged.txt examples/Berlin.ectopics.tsv`

User can run this test by executing *runtest.sh* script:
`./runtest.sh`

