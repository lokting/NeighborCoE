# NeighborCoE

NeighborCoE is a Python-based pipeline for identifying operons from bacterial genome annotation and count data, and for further defining the core operonome and pan operonome.

## Features
- Parse GBFF/GFF annotations into operon-ready tables
- Integrate count files and genomic context
- Infer operons based on intergenic organization
- Annotate operon gene composition for each strain
- Identify the core operonome across strains using cluster-mapped operon signatures

# PanOperonome
Automated tool for pan-operonome network analysis.
  
## When using this package
## Some input files are necessary

1.Count File

2.Gff Gbff File

3.Chromosome and plasmid ID (including its strand)

4.Genome length

5.Core-genome cluster file for core operonome identification

6.Pan-genome results and mutual blast results for pan operonome identification

The example file format is in images/
