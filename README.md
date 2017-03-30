# Pedselect v0.1
Pedselect is a Python script for selecting and prioritizing individuals to sequence from large multiplex pedigrees

## Overview
Large multiplex pedigrees can be informative for helping identify disease associated risk loci. As many pedigrees are collected over long periods of time, pedigree members may vary in terms of whether their phenotype has been observed, whether they have already been sequenced/genotyped, and whether they have DNA collected that can readily be sequenced/genotyped. Given budget constraints, only sequencing some members of the pedigree may be feasible. Pedigree-based genotype imputation methods (e.g. Cheung et al., 2013) can then be used to fill-in the genotypes of members who have not (or cannot) been sequenced.

Pedselect is a tool for selecting the most informative individuals to perform whole-genome-sequencing/genotyping in a pedigree. In order to maximize power to identify disease associated loci, Pedselect first prioritizes members whose phenotypes have been collected and have DNA available. It then sequentially selects members to sequence based on their ability to inform the genotypes of nearby relatives. Once all individuals with phenotypes and available DNA have been selected, Pedselect then sequentially selects among the remaining individuals who have DNA available but have not been phenotyped based on their ability to inform the genotypes of remaining phenotyped individuals (who do not have DNA available).

For an overview of the method, see Pedselect.pdf.

## Usage
Pedselect requires Python 2.7. To clone respository:
```
git clone https://github.com/nygenome/pedselect
```

### Input file
A space-delimited file with one individual per-line and the following columns:
```
Individual ID
Father's ID
Mother's ID
Sex
Whether individual has already been sequenced/genotyped (1 or 0)
Whether individual can be sequenced sequenced/genotyped (1 or 0)
Whether individual has been phenotyped (1 or 0)
```

To run, type:
```
python pedselect.py -p example.ped -m 6
```

The optional -m flag denotes the number of meiosis events from an individual to consider when calculating that individual's imputability score.

### Output
The script will print out a list of individuals ordered by priority for sequencing. Next to each individual is an overall imputability score for the rest of the pedigree assuming that individual and all those preceeding them has been sequenced.

## License
This project is licensed under BSD License 2.0

## Contact
Jimmy Liu (firstname dot z dot lastname dot gmail dot com), New York Genome Center

March 30 2017


