
# Homology segment analysis

The application estimates the DNA sequence homology, finding short 
fragments (k-mers) that match for the two extended genomic areas. 
It also calculates correlations between 
the frequencies of matched fragments (FMF) and any other genome data. 
In example data the frequencies of ectopic contacts (FEC) between the 
Drosophyla polytene chromosome sections are used.

## The software capabilities

* Finds k-mers of the provided DNA sequence with genome and computes their 
frequencies (frequencies of matched fragments FMF).
* Calculates a correlation between FMF and numeric data specified for genome (for example, FEC).
* Compare the correlation levels and proportions of significant correlations calculated at different paramaters.
* Plots homology level graph.

## Licence
This program is licensed under GNU GPLv3 licence: https://www.gnu.org/licenses/gpl.html

# Installation 
The software is written using Python language version 3.
In addition to main Python package, it uses the following Python modules:

* Aho-Corasick for fast fragments search.
* Numpy and Scipy for correlation analysis.
* Matplotlib for plotting graphs.

Both for Windows and Linux users we recommend installing all dependencies via Miniconda. 
The same way should work for Mac also. 
However, the application was not tested on Mac.

## Installing Anaconda
* Download installer for the latest Python3 version from miniconda site: 
https://conda.io/en/latest/miniconda.html 
* Follow the installation instructions.

## Creating Conda environment
Execute the following commands from your shell

    conda create --name py3 python=3.8
    conda activate py3
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda update --all
    conda install numpy scipy matplotlib pandas pyahocorasick

### Windows users note
If installation of pyahocorasick module fails, try to install 
Visual Studio build tools from here: 
https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017
Note: currently, there may be problems with using analysis.py script with Windows. 

### Linux users note
All Linux distributions contains Python3 by default. 
So Linux users instead of Anaconda package can install all required packages 
via distribution package manager and python pip.

Example for Ubuntu/Debian:
    sudo apt install build-essential python3-dev
    pip3 install numpy scipy matplotlib pandas pyahocorasick

## The main program installation
Just checkout from master branch.

`git clone https://bitbucket.org/beneor/homology-segment-analysis.git`

Then navigate to program base folder and start.

# Usage

## Calculation of FMF for provided DNA sequence (segment) and the whole 
genome using segmentanalysis.py

The `segmentanalysis.py` script splits DNA segment into a set of short 
(about 10 - 60 nt) fragments (k-mers), then splits the whole chromosome into a set 
of long intervals (chunks; 10 kb by default), and seeks for k-mers in 
the chunks. Then it computes the frequencies 
of match normalized by chunks (the average FMF is 1 for all chunks) and 
calculates k-mers matching frequencies for the long chromosome areas, 
if the information about areas coordinates is provided. 
By default: areas are the Drosophila polytene chromosome cytoband, 
or sections (X1A - 4102E; file containing sections borders: 
`./examples/DmelMapTable.160615c.bed`) 
or zones without division into sections (X1 - 4102; file containing zone 
borders: `./examples/DmelMapTable.2.bed`).
Resulting data are grouped by chunks and by cytobands.

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
File name can be omitted to take the segment from the same genome as the analyzed one.
Start and stop locations also can be omitted to use the whole FASTA record. 

##### Examples
* `examples/dm6.onlyX.fa.gz :X:103614:408582` -- part of Drosophila melanogaster X chromosome sequence from file examples/dm6.onlyX.fa.gz (positions 103614:408582)
* `examples/dm6.onlyX.fa.gz` -- the whole Drosophila melanogaster X chromosome sequence sequence from file examples/dm6.onlyX.fa.gz

### Optional arguments
* `-h, --help`  -          Show help message and exit
* `-v, --verbose`  -       Print additional information to stdout
* `-s, --fragmentsizes` - Set of fragment sizes to search divided by comma. Default: 20,25,30,35,40
* `-d, --fragmentdensity`  - Use random chosen fragments instead of all possible.
                        Sets average frequency of fragments in letters. For
                        example 5 means that each 5-th letter will be start of
                        fragment
* `-c, --chunk`  - Chunk size to divide chromosome, in kilobases. Default: 10. 
                       
* `--nodump` - Do not save found fragments
* `--mindumpsize` - Minimum size of fragment to store exact locations in file. Default 20
* `-b, --blacklist` - List of fragment to forcedly exclude from search. 
Can be used to exclude overpresented repeats. The example file containing repeats: `./repeats.txt`.
* `-i, --include` - List of fragments to forcedly include in search.  
The example file containing fragments to include is: `./include.txt`. Note that the excluede and included sequences must be 
of the same length as the length of fragments. 

### Output data

Folders containing files with FMF data and fragments are written to the
same folder, where input genome file is located. The folder name is constructed from genome file name and segment.
For example, `dm6.onlyX.fa.gz.X-103614-408582` folder containes FMF 
data and fragments for 103614-408582 bp area of the _Drosophila_ X-chromosome. 

#### Output files naming

* `cytocounts.lN-A-B.txt` - files contining normalized FMF values for 
chunks united to sections (cytobands).
* `ncounts.lN-A-B.txt` - files containing normalized FMF valus for chunks
 not united to cytobands.
* `fragments.lN-A-B.txt`- files containing matching fragment sequences and 
their locations in genome.

* `lN` - fragment length (e.g. `l30` - 30 nt).
* `A` - the cytoband start position (e.g. `103614`)
* `B` - `dir`, `rev` or `merged` (fragments of direct, reverse or both DNA 
strand of the chromosome segment, respectively)

The cytocount and ncount files content: chromosome name, start, stop, 
cytoband or chunk name, normalized FMF.
The fragment files content: chromosome name, fragment sequence, locations. 
The content is arranged by the fragment sequnces.

### Usage example

`./segmentanalysis.py -s 30,50 -b repeats.txt -i include.txt examples/dm6.onlyX.fa.gz :X:11982050:12772075 examples/DmelMapTable.160615c.bed`
Run FMF calculation for _Drosophila_ X-chromosome 11AB region

## Data preparation for analysis of FMF-FEC correlations.
The pre-calculated matrix of FMF must be built before running `analysis.py` script.

### Calculating FMF for whole chromosome
The `fmf_chromosome.sh` script is used to calculate FMF for one particular chromosome each chromosome bands.
The script running format is following:

`./fmf_chromosome.sh <genome> <cytobands> <chromosome>`

It creatsd the batch file to process all bands and runs it in parallel the `segmentanalysis.py` script.
Other run parameters (blacklist, fragment sizes) are hardcoded in the header section of the script.

#### Usage example

 `./fmf_chromosome.sh examples/dm6.onlyX.fa.gz examples/DmelMapTable.160615c.bed X`
 
 Calculation for _Drosophila_ X chromosome for two fragment lengths takes about 5-10 minutes on 6-cores CPU.

Run parallel FMF calculations for Drosophila X-chromosome
(for all X chromosome disc=sections, 
no fragments excluding, fragment lengths 30 and 50 nt).

### FMF data processing

The script `Filesmodification.py` in `./examples` folder must be run to 
arrange FMF and fragment files for further analysis,
placing them to appropriate folders of the ./cytocountData folder 

Running the script from the `./examples` folder: `python3 Filesmodification.py`

The user is asked to specify the type of the folder depending on the type 
of computations performed at the previous stage: 
`discs` - with division into sections; 
`zones` - without division into sections; 
`all` - all fragments, `nr` - no repeats. 
 
The FMF and fragment files are placed to the folder named 
according to the user choice:
Example: `discs/all/l30-cytocounts` - the folder containing FMF for all 
matching fragments of 30 nt length calculated for the chromosome sections; 
`discs/all/l30-fragments` -  the corresponding fragment sequences.

The resulting folders is written (or must be manually placed) to 
`./example/Cytocountsdata` folder. The folder should have the following structure:
    
    Cytocountsdata
        - discs
            - all
                -l<length>-cytocounts
                -l<length>-fragments
            - nr
                -l<length>-cytocounts
                -l<length>-fragments
        - zones
            - all
                -l<length>-cytocounts
                -l<length>-fragments
            - nr
                -l<length>-cytocounts
                -l<length>-fragments


The file names in the folder consists from the zone number and disc=section number.
For the _Drosophila melanogaster_ X chromosome, there are 19 zones and 6 sections (A-F) in each zone. 
Note that the section notation changes in the new file names: 
"A, B, C, D, E, F" are replaced by "1, 2, 3, 4, 5, 6", so "1A" is "11" now, "13F" is "136", etc. 
For example, the file 123.txt corresponds to the section 12C.

The example FMF and fragments data precalculated for the whole X chromosome 
(sections=discs, all fragments, fragment lengths 30 and 50) 
are given in `Cytocountdata-examples.tar.gz` archive.

## Ectopics data (FEC)
Currently FEC are organized in set of 1-dimensional file tables.
The data format change to more convenient 2D-matrix is under development. 

The ectopic data folder should have the following structure:

    Xectopicsdata
        - strain 1
            - discs
            - zones
        - strain 2
            - discs
            - zones

The file names follows the same structure as for cytocounts FMF data
* The file name for zones is just a number of zone
* The file name for section consists of the zone number (1-19) and the 
number of section (1-6). So the file 123.csv represents the section 12C

The sample data for _Drosophila melanogaster_ X-chromosome are given in 
`Xectopicsdata` folder, being also available here: 
https://drive.google.com/file/d/1cPDRCaTrUlxa4Pp0fR4DxnP9NzEKLCgt/view?usp= 

## Full-chromosome correlation analysis using analysis.py 

The script `analysis.py` computes the average values of correlations 
between FMF and any other genome data (by default: FEC), comparing them 
for different strains and fragment length:

1.) Spearman rho specific and unspecific correlation values are calculated 
using FEC and FMF data obtained for the whole X chromosome. 
Specific correlation is the correlation between FEC and FMF data obtained 
for the same pair of sections, unspecific correlation - for the different pairs of sections.

2.) For statistically significant correlations: calculation of the average 
rho values (R) and the proportion of significant correlations (P). 
The data are saved in: 
`./Correlation/Group_name/Group_name-60-analysis-average-correlations.xls`. 
The files content: the average rho (R) and P values and their statistical 
errors, calculated at all fragment lengths. 
Group_name is notated as A-B-C-D; A - strain, B - discs or zones, 
C - all fragments (all) or no repeats (nr), D -  specific or unspecific correlation. 
All possible combinations of parameters are automatically varied. 

3.) Comparison of R values using the Mann-Whitney U-test and P values 
using the proportion test. 
R and P are compared for different fragment lengths or for the same 
fragment length with different parameter sets (Group_name). Thereby all 
the compared groups differ in exactly one parameter.
The data are saved in `./MW-analysis/diff-length/Correlations-Group_name1-Group_name2-MW-test.xls` 
and `./MW-analysis/diff-parameters/Correlations-Group_name1-Group_name2-MW-test.xls`, respectively.
The files content: name 1 and name 2 - the fragment lengths for Group 1 
and Group 2, prob_cor and prob_p - statistical sinificance of hypothesis 
that the values are the same for R values and P values, respectively. 
Only the data for prob_cor or prob_p < 0.05 are given. R and P distributions 
are assayed as well: norm - normal distribution , not-norm - the distribution is not normal. 

4.) For FEC and FMF data obtained for all the X chromosome parts of a 
specific length (D): calculation of the average R and P values. 
Calculations are performed for specific paramenter sets (Group_name). 
The data are saved in: `./Nearest/Group_name/l-nearest-av-y.xls`, 
l - the short fragment length, y - specific or unspecific correlation.
The files content: area - D length (in sections), average_correlation 
and average share - R and P values averaged for all the chromosome 
parts of the D length.

5.) For two different groups (1 and 2) and D: comparison of the average 
R values using the Mann-Whitney U-test and the average P values using 
the proportion test.
The data are saved in: `./MW-analysis/nearest/Archive/Nearest-Group_name1-Group_name2-l-avcorr.xls`; 
1 - Group_name 1, 2 - Group_name 2, l - the short fragment length.
The files content: segmentlength - D (in sections), avcor_p and 
avprob_p - R and P values averaged for all chromosome segments of 
the D length, corst_err1, corst_err2 - statisticall errors of 
the averaged R and P for groups 1 and 2. Only the data for 
prob_cor or prob_p < 0.05 are given. R and P distributions are assayed as well.
 
General run syntax is following: 
`./analysis.py <Cytocountsdata> <Xectopicsdata> [--nearest]` 

The main arguments are:
+ `Cytocountsdata` - Folder containing pre-calculated data with cytocounts.
+ `Xectopicsdata` - optional tab-delimeted file containing ectopic contacts data

### Optional arguments
+ `strain1` - Name of the first strain to analyse. Default: _CS_
+ `strain2` - Name of the second strain to analyse. Default: _agn_
+ `-h, --help` -     Show help message and exit
+ `-s, --fragmentsizes` -- Set of fragment sizes to search. Default: 30,50
+ `-N,  --nearest` - Perform nearest calculations. Default: No The part 
of the application is currently not refactored, so using it is not so fast as other scripts.

### Output data examples
After the analysis is completed, the data folder should have the following structure:

    - Correlations
        - folder Group_name (e.g. 'CS-discs-all-specific': specific correlation for CS strain, discs=sections, all fragments)
            - Discs
                -working files *.csv3 containing only statistically significant Spearman rho values and significaтсу levels, for fragments of specific length (e.g. l10: 10-specifcor.csv3)
            - Archive
                -working files *.csv2 containing Spearman rho values and significaтсу levels for all sections, for fragments of specific length (e.g. l10: 10-specifcor.csv2)
            - files Group_name-60-analysis-average-correlations.xls
    - MW-analysis
        - diff-length
            - files Correlations-Group_name1-Group_name2-MW-test.xls
        - diff-parameters
            - files Correlations-Group_name1-Group_name2-MW-test.xls
        - nearest
            - Archive
                    - files Nearest-Group_name1-Group_name2-l-avcorr.xls
    - Nearest
        - folder Group_name (e.g. 'CS-discs-all-specific': specific correlation for CS strain, discs=sections, all fragments)
            - files *.l-nearest-av-y.xls

The example of the analysis data is given in `./examples/analysis-examples.tar.gz` archive.           

### Run example
`./analysis.py examples/Cytocountsdata examples/Xectopicsdata/`


## Generation of arranged fragment list
The script `Fragments-counts.py` in `./Additional scripts` folder should be used to generate the list of short matching fragment 
of a specific length arranged according to the number of their occurrence in genome. 
Place this script in the folder containing fragment sequences (`-l<length>-fragments`) and run there: 
`python3 Fragments-counts.py`
The result file `Result.csv` can be opened and analyzed using LibreOffice Calc. The file content: fragent number of occurrence (>10), fragment sequence.

## Homology graph plotting using genomecorrelation.py 

The `genomecorrelation.py` script plots homology level graphs (both grouped by chunks or by cytobands),
and calculated correlation of homology level with experimental frequencies, if provided (by default: FEC).

General run syntax is the following: 

`genomecorrelation.py <options> homologyFileName [ectopicsFileName]` 

The main arguments are:

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


## Test run examples
Generate homology data for Drosophila genome, included in examples data set 
`./segmentanalysis.py -v -s 30 -d 2.5 examples/dm6.nounmapped.fa.gz :X:11982050:12772070 examples/DmelMapTable.160615c.bed`

View homology graph together with ectopics contacts frequency
`/genomecorrelation.py -m -c X,2L examples/dm6.nounmapped.fa.gz.X-11982050-12772070/cytocounts.l30-merged.txt examples/Berlin.ectopics.tsv`

User can run this test by executing *runtest.sh* script:
`./runtest.sh`

