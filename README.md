# Hi-Cluster
Hi-Cluster is a method that generate multiple MAGs (metagenome-assembled genomes) using Hi-C reads.

# Usage

Hi-Cluster is composed of two scripts (link_counter.py and Hi-Cluster.py) and can be devided into two steps.

## step1 (link_counter.py)
```
python link_counter.py samfile.sam scaffolds.fa cutting_site outputFileName minimum_length scaffoldList.txt

samfile.sam       Mapping result of Hi-C reads to scaffolds. The file should be sam format.
scaffolds.fa      input fasta file of scaffolds.
cutting_site      restriction site (e.g. MluCI -> AATT).
outputFileName    out put file.
minimum_length    minimum scaffold length.
scaffoldList.txt  header of samfile. (which contain only a scaffold name per line.)
```
