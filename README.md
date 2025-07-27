# Double Gyrase Analysis Pipeline

This repository contains scripts and data for analyzing double gyrase gene sequences and their mutations.

## Project Structure
```
double_gyrae/
├── input/             # Input FASTA files and logits
├── figures/           # Generated plots and visualizations
├── scripts/           # Analysis scripts
├── output/            # Generated output files
├── run_evo.sh        # EVO model execution script
└── README.md
```

## Data Generation Pipeline

### 1. Sequence Preparation
- Initial sequences were prepared in FASTA format:
  - `double_ref.fasta`: Reference double gyrase sequence
  - `double_83_TTT`: Mutant sequence with TTT mutation
  - `double_random`: Random sequence of same length

### 2. EVO Model Analysis
From the project root directory, run:
```bash
# Run the EVO model
./run_evo.sh double-genes input/double_genes.fasta input
```
This will:
- Create necessary directories (output, jobs, figures)
- Submit the job to the EVO model
- Download results when complete
- Generate logits files in the format: `input/input_[sequence_name]_logits.npy`

### 3. Data Analysis and Visualization
The analysis is performed using R scripts in the following order:

1. **Initialize Data Processing**
```R
scripts/handle_EVO2_output_v1.R  # Handles logits file processing
scripts/evo2_analysis_functions_v1.R  # Core analysis functions
```

2. **Main Analysis Script**
```R
scripts/EVO_learned_resist.R  # Generates all figures
```

This script performs:
- Loads and processes logits data
- Calculates log-likelihoods
- Generates difference plots between sequences
- Creates visualizations for:
  - Reference vs Random comparisons
  - Mutant vs Random comparisons
  - Second gene region analysis
  - Transition region analysis

## Generated Figures

The analysis produces several PDF figures in the `figures/` directory:
- `double_ref-random.pdf`: Comparison of reference to random sequence
- `double_mut-random.pdf`: Comparison of mutant to random sequence
- `83_TTT_2nd.pdf`, `83_ref_2nd.pdf`: Analysis of second gene region
- `trans_*.pdf`: Analysis of transition regions

## Dependencies

- R packages:
  - Biostrings
  - ggplot2
  - reticulate (for Python integration)
- Python with numpy

## Usage

1. Clone the repository
2. Place your FASTA files in the `input/` directory
3. From the project root directory, run the EVO model:
```bash
./run_evo.sh double-genes input/double_genes.fasta input
```
4. Run the R analysis script:
```R
source("scripts/EVO_learned_resist.R")
```

## Notes

- This analysis of double genes indicates that evo2 7b memorizes the gene 
sequence within ~25bp and can then accurately predict every bp of the repeated
gene going forward
- The highlighted nucleotides are the frist 3 nucleotides of the repeated gene
