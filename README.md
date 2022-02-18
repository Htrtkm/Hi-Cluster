# Hi-Cluster

Hi-Cluster is an method that generate multiple MAGs (metagenome-assembled genomes) using Hi-C reads.

# Usage

Hi-Cluster is composed of three scripts (link_counter.py, merge_contactfiles.py and Hi-Cluster.py) and can be devided into three steps.

## step1 (link_counter.py)
```
python link_counter.py (1)samfile.sam (2)scaffolds.fa (3)cutting_site (4)outputFileName (5)minimum_length (6)scaffoldList.txt

samfile.sam       Mapping result of Hi-C reads to scaffolds. The file should be sam format.
scaffolds.fa      input fasta file of scaffolds.
cutting_site      restriction site (e.g. MluCI -> AATT).
outputFileName    out put file.
minimum_length    minimum scaffold length.
scaffoldList.txt  header of samfile. (which contain only a scaffold name per line.)
```

## step2 (merge_contactfiles.py)
