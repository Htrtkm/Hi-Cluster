# Hi-Cluster

Hi-Cluster is an method that generate multiple MAGs (metagenome-assembled genomes) using Hi-C reads.

# Usage

Hi-Cluster is composed of three scripts (link_counter.py, merge_contactfiles.py and Hi-Cluster.py) and can be devided into three steps.

## step1 (link_counter.py)
this script generates contact information from samfiles.
```
python link_counter.py (1)samfile.sam (2)scaffolds.fa (3)cutting_site (4)outputFileName (5)minimum_length (6)scaffoldList.txt

samfile.sam       Mapping result of Hi-C reads to scaffolds. The file should be sam format.
scaffolds.fa      input fasta file of scaffolds.
cutting_site      restriction site (e.g. MluCI -> AATT).
outputFileName    output file.
minimum_length    minimum scaffold length.
scaffoldList.txt  header of samfile. (which contain only a scaffold name per line.)
```

## step2 (merge_contactfiles.py)
this script integrates two contact files and generate a gml file which is information of proximity among scaffolds with the style of network.
```
python merge_contactfiles.py (1)contactFile1.txt (2)contactFile2.txt (3)outputGMLfile.gml (4)scaffoldList.txt (5) minimum_length

contactFile1      output of link_counter.py
contactFile2      output of link_counter.py
outputGMLfile     output of merge_contactfiles.py
scaffoldList.txt  header of samfile. (which contain only a scaffold name per line.)
minimum_length    minimum scaffold length.
```

## step3 (Hi-Cluster.py)
this script creates MAGs from gml file.
```
python Hi-Cluster.py (1)inputGMLfile.gml (2)fastafile.fa (3)threashold (4)outputFileName

inputGMLfile.gml  input gml file
fastafile.fa      input fasta file
threashold        default (0.96)
outputFileName    output file name
```
