
# Homology Modeling For 3D Protein Structure Prediction

Welcome to the Protein Structure Prediction GitHub repository! Our project focuses on advancing the understanding of protein structures by employing computational methods, specifically comparative modeling. We recognize the significance of determining 3D protein structures for applications ranging from drug design to unraveling disease mechanisms. This toolkit utilizes the RCSB Protein Data Bank, incorporating AlphaFold, to retrieve amino acid sequences. Leveraging the Basic Local Alignment Search Tool (BLAST) and efficient sequence alignment techniques, we enhance the accuracy of predictions. Our post-processing methods optimize alignments, and PDB file generation produces organized files for further analysis. The project culminates in a comparative analysis with AlphaFold, ensuring the reliability of our predictions. 

## Getting Started

### Requirements

Before you begin, make sure you have Python installed on your system. You can download Python from [python.org](https://www.python.org/).

### Installing Dependencies

To run this project, you'll need to install the following Python modules. You can do this using the following `pip` commands:

```bash
pip install requests
pip install numpy  # This may be required by Biopython
pip install biopython
pip install pard  # Make sure to check for the latest version on PyPI
```
## Usage

### BLAST

BLAST (Basic Local Alignment Search Tool) is used to compare protein query sequences with a protein database to identify homologous proteins.

To run BLAST:

```bash
python blast_new.py
```