
# Homology segment analysis

These scripts analyse homology level of particular sequence with genome, plot homology level graph 
and calculated correlation between homology level any some genome-related data (initially level of ectopic contacts). 
Resulting data are grouped by genome intervals (chunks) and by cytobands (if cytomband data is avaliable)

First script `segmentanalysis.py` splits segment into set of short (20 - 40 kb) fragments,
then splits chromosome into a set of long (10 kb by default) chunks, and seeks for positions where these fragments matches with genome.
Then it computes frequencies of occurense for segments within chromosome,
makes table of the match frequencies normalized by chunks, and calculates fragments matching frequency by cytobands, 
if cytoband information is provided.

Second script `genomecorrelation.py` plots homology level graphs (both grouped by chunks or by cytobands),
and calculated correlation of homology level with experimental data, if provided.


# Prerequirements
Software is written using Python language version 3.
In addition to main Python package, it uses following Python modules:
* Aho-Corasick for fast fragment search
* Numpy and Scipy for correlation analysis
* Matplotlib for plotting graphs.

## Installation for windows  
For Windows users we recommend installing Anaconda package, containing all required
Python packages. Visit https://www.anaconda.com/ and follow installation instructions.

Next, install visual studio buildtools form here: https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017

Finally open Anaconda console and install PyAhoCorasick library

`pip install pyahocorasick`

## For Linux

Install Anaconda package, containing all required
Python modules. Visit https://www.anaconda.com/ and follow installation instructions.

Next, install base development packages via distribution package tool.
Please refer to distribution documentation for it.
For example for Ubuntu command should be like this:

`sudo apt install build-essential python3-dev`

Finally install PyAhoCorasick library using pip

`pip3 install pyahocorasick`

#### Note 
All Linux distibutions contains Python3 by default. 
So Linux users instead of Anaconda package can install all required packages using following command:

`pip3 install numpy scipy matplotlib`

## Main program installation

Just checkout from master branch.

`git clone https://bitbucket.org/beneor/homology-segment-analysis.git`

Then navigate to program base folder, and start using it.

# Usage

## Homology calculation using segmentanalysis.py 

### Main arguments
General run syntax is following: 

`segmentanalysis.py <options> chromosomeseq segment [cytobands]`

The main arguments are

+  `fastaFileName` -      FASTA file (may be gzipped) containing genome sequence                      
+  `segment` -              segment of chromosome to analyze: file and location. See format notes below.
+  `cytobands` - optional BED-format file containing locations of cytobands.
Or any other intervals to group homology data 

#### Segment notation
**segment** option has flexible format to handle multiple cases. 
In general it contains name of FASTA file and location of segment in this file in BED 
notation (chromosome, start,stop) divided by **:**
File name can be omitted to take segment from same genome as analyzed one.
Start and stop locations also can be omitted to use whole fasta record. 

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

## Test run examples
Generate homology data for Drosophila genome, included in examples data set 

`./segmentanalysis.py -v -s 30 -d 2.5 examples/dm6.nounmapped.fa.gz :X:11982050:12772070 examples/DmelMapTable.160615c.bed`

View homology graph together with ectopics contacts frequency

`/genomecorrelation.py -m -c X,2L examples/dm6.nounmapped.fa.gz.X-11982050-12772070/cytocounts.l30-merged.txt examples/Berlin.ectopics.tsv`

User can run this test by executing *runtest.sh* script:

`./runtest.sh`

## Licence
This program is licensed under GNU GPLv3 licence
https://www.gnu.org/licenses/gpl.html 

</main>
