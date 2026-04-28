# NeighborCoE

NeighborCoE is a Python-based toolkit for bacterial operonome analysis. It includes three main modules:

1. **GetOperon**: identification of operons from genome annotation and transcriptomic count data.
2. **CoreTUs**: identification of core transcriptional units / core operonome across strains.
3. **PanTUs**: construction of pan-operonome networks from strain-specific operon annotations and homologous-gene relationships.

Representative input and output files are provided in the `examples/` directory.

---

## Repository structure

```text
NeiborCoE/
├── GetOperon/
│   ├── code/
│   │   ├── 0.gbff_to_excel.py
│   │   ├── 1.gff_to_excel.py
│   │   ├── 2.count_file.py
│   │   ├── 3.get_intergenegap.py
│   │   ├── 4-0.get_operon.py
│   │   └── 4-1.get_operon_function.py
│   └── operon.sh
├── CoreTUs/
│   └── get_core_operonome.py
├── PanTUs/
│   ├── main.py
│   └── src/
│       ├── stage_0_input.py
│       ├── stage_1_preprocess.py
│       ├── stage_2_popgid.py
│       ├── stage_3_network.py
│       └── stage_4_finalize.py
├── examples/
│   ├── GetOperon/
│   │   ├── GetOperon_input_data/
│   │   └── GetOperon_output_data/
│   ├── CoreTUs/
│   │   ├── CoreTUs_input_data/
│   │   └── CoreTUs_output_data/
│   └── PanTUs/
│       ├── PanTUs_input_data/
│       └── PanTUs_output_data/
├── requirement.txt
└── README.md
```

---

## Installation

Clone the repository and enter the working directory:

```bash
git clone https://github.com/lokting/NeighborCoE.git
cd NeighborCoE
```

Create a conda environment:

```bash
conda create -n NeighborCoE python=3.8 -y
conda activate NeighborCoE
```

Install the required Python packages:

```bash
pip install -r requirement.txt
```

---

## Module 1. GetOperon

### Purpose

`GetOperon` identifies operons using genome annotation files, read count files, intergenic distance, gene orientation, and expression correlation between adjacent genes.

### Input files

Place all input files under:

```text
GetOperon/input/
```

The required input structure is:

```text
GetOperon/input/
├── Count_file/
│   ├── sample1_count.txt
│   ├── sample2_count.txt
│   └── ...
├── genome_gbff/
│   └── genome.gbff
├── genome_gff/
│   └── genome.gff
├── chr_list.txt
├── pla_list.txt
├── genome_length.txt
└── local_tag_dele.txt
```

Example files are provided in:

```text
examples/GetOperon/GetOperon_input_data/
```

### Run GetOperon

Enter the `GetOperon` directory:

```bash
cd GetOperon
```

Prepare the input and output folders:

```bash
mkdir -p input output tmp_output
```

If you want to run the provided example data, copy the example input files into `input/`:

```bash
cp -r ../examples/GetOperon/GetOperon_input_data/* input/
```

Run the full GetOperon pipeline:

```bash
bash operon.sh
```

Return to the main repository directory:

```bash
cd ..
```

### Output files

The main output file is:

```text
GetOperon/output/operon-gene.xlsx
```

This file contains the predicted operon ID, gene composition, and functional annotation.

---

## Module 2. CoreTUs

### Purpose

`CoreTUs` identifies core operons / core transcriptional units shared across multiple strains based on core-genome clusters and strain-specific operon annotations.

### Input files

Place all input files under:

```text
CoreTUs/input/
```

The required input structure is:

```text
CoreTUs/input/
├── CG_ALL.txt
└── operon_file/
    ├── strain1_operon-gene.xlsx
    ├── strain2_operon-gene.xlsx
    └── ...
```

Example files are provided in:

```text
examples/CoreTUs/CoreTUs_input_data/
```

### Run CoreTUs

Enter the `CoreTUs` directory:

```bash
cd CoreTUs
```

Prepare the input and output folders:

```bash
mkdir -p input output
```

If you want to run the provided example data, copy the example input files into `input/`:

```bash
cp -r ../examples/CoreTUs/CoreTUs_input_data/* input/
```

Run the CoreTUs script:

```bash
python get_core_operonome.py
```

Return to the main repository directory:

```bash
cd ..
```

### Output files

The main output file is:

```text
CoreTUs/output/core_operonome.txt
```

This file contains the identified core operonome across the input strains.

---

## Module 3. PanTUs

### Purpose

`PanTUs` constructs pan-operonome networks using strain-specific operon annotations and homologous-gene / pan-genome information.

### Input files

Place strain-specific operon files under:

```text
PanTUs/data/input/
```

Place reference files under:

```text
PanTUs/data/reference/
```

The required input structure is:

```text
PanTUs/data/
├── input/
│   ├── strain1_operon-gene.xlsx
│   ├── strain2_operon-gene.xlsx
│   └── ...
└── reference/
    ├── PG.txt
    └── nPO_homo.xlsx
```

Example files are provided in:

```text
examples/PanTUs/PanTUs_input_data/
```

### Run PanTUs

Enter the `PanTUs` directory:

```bash
cd PanTUs
```

Prepare the input, reference, and output folders:

```bash
mkdir -p data/input data/reference output
```

If you want to run the provided example data, copy the example input files into `data/`:

```bash
cp -r ../examples/PanTUs/PanTUs_input_data/input/* data/input/
cp -r ../examples/PanTUs/PanTUs_input_data/reference/* data/reference/
```

Run the full PanTUs pipeline:

```bash
python main.py
```

Return to the main repository directory:

```bash
cd ..
```

### Output files

The main output file is:

```text
PanTUs/output/network_results/network_final.csv
```

This file contains the final pan-operonome network.

---

## Notes

- All scripts should be run from their own module directory. For example, run `bash operon.sh` inside `GetOperon/`, not from the main `NeighborCoE/` directory.
- The example input and output files are provided only as representative data to show the required file formats.
- The required input files must be kept in the input folder even if they do not contain any records. Do not delete empty input files; keep them as blank placeholder files when no corresponding data are available.

