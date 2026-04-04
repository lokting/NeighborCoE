# Operonome-NeighborCoE

NeighborCoE is a Python-based toolkit for bacterial operon analysis, including:
1. operon identification from genome annotation and count data,
2. core operonome identification across strains, and
3. pan-operonome network construction.

## Overview

This repository currently contains two modules:

- **Operon/Core Operonome module**: identifies operons from genome annotation and count data, and further defines the core operonome across strains.
- **PanOperonome module**: constructs pan-operonome networks from strain-level operon annotations and pan-genome / homologous-gene relationships.

## Repository structure

```text
Operonome/
├── get_operon/
│   └── scr/
│       ├── 0.gbff_to_excel.py
│       ├── 1.gff_to_excel.py
│       ├── 2.count_file.py
│       ├── 3.get_intergenegap.py
│       ├── 4-0.get_operon.py
│       ├── 4-1.get_operon_function.py
│       └── 5.get_core_operonome.py
├── PanOperonome/
│   └── src/
│       ├── stage_0_input.py
│       ├── stage_1_preprocess.py
│       ├── stage_2_popgid.py
│       └── stage_3_network.py
├── images/
├── README.md
└── requirement.txt

## Module 1. Operon and core operonome

### Input files

Required inputs include:

- Count file
- GFF/GBFF annotation file
- Chromosome and plasmid IDs (including strand information)
- Genome length
- Core-genome cluster file for core operonome identification

Example input formats are provided in `images/`.

### Workflow

Run scripts in the following order:

1. `0.gbff_to_excel.py`
2. `1.gff_to_excel.py`
3. `2.count_file.py`
4. `3.get_intergenegap.py`
5. `4-0.get_operon.py`
6. `4-1.get_operon_function.py`
7. `5.get_core_operonome.py`

### Output

Typical outputs include:

- processed annotation tables
- intergenic-gap tables
- strain-specific operon annotation tables
- core operonome results

## Module 2. PanOperonome

### Input files

Required inputs include:

- strain-specific operon-gene annotation files
- pan-genome cluster / PG mapping file
- homologous operon summary or mutual BLAST-derived operon relationship file

### Workflow

Run scripts in the following order:

1. `stage_0_input.py`
2. `stage_1_preprocess.py`
3. `stage_2_popgid.py`
4. `stage_3_network.py`

### Output

Typical outputs include:

- pan-operonome network results

## Notes

- Scripts are designed to be run sequentially.
- Please adjust input/output paths according to your local environment.
- Example file formats are provided in `images/`.
