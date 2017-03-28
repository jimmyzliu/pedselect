# Pedselect
Pedselect is a Python script for selecting and prioritizing individuals to sequence from large multiplex pedigrees

## Overview
Large multiplex pedigrees can be informative for indentifying disease associated risk loci. As many pedigrees are collected over long periods of time, pedigree members may vary in terms of whether their phenotype has been observed, whether they have already been sequenced/genotyped, and whether they have DNA collected that can readily been sequenced/genotyped. Given budget constraints, only sequencing some members may be feasible. Pedigree-based genotype imputation methods can then be used to 'fill-in' the genotypes of members who have not (or cannot) been sequenced.

Pedselect is a tool for selecting the most informative individuals to perform whole-genome-sequencing/genotyping in a pedigree. In order to maximize power to identify disease associated loci, Pedselect first prioritizes members whose phenotypes have been collected and have DNA available. It then sequentially selects members to sequence based on their ability to inform the genotypes of nearby relatives. Once all individuals with phenotypes and available DNA have been selected, Pedselect then sequentially selects among the remaining individuals who have DNA available (but have not been phenotyped) based on their ability to inform the genotypes of remaining phenotyped individuals (but who do not have DNA available).

## Usage
Pedselect requires Python 2.7. To clone respository:
```
git clone https://github.com/jimmyzliu/pedselect
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

The optional -m flag denotes the number of meiosis events from an individual to consider when calculating that individual's imputability score (see Method).

### Output
The script will print out a list of individuals ordered by priority for sequencing. Next to each individual is an overall imputability score for the rest of the pedigree assuming that individual and all those preceeding them has been sequenced (see Method).

## Method
In order to quantify “imputability”, we derived a score for each individual i for whom we wish to impute: , where j = 1,2,…n are the genotyped relatives of individual i and m is the number of meiosis events between individual i and j. When calculating Si, we do not count a genotyped distant relative if a closer genotyped relative along the same relatedness path is observed. For example, if the genotypes of both the mother and maternal grandmother are available, we do not consider information from the grandmother since all shared identical-by-descent segments with the grandmother are also observed in the mother (Rafnar et al., 2011).

We use this imputability score to prioritize which individuals with available DNA and a known phenotype to sequence. Our algorithm goes through each of these individuals in turn and, assuming that the individual is sequenced, calculates the imputability score of the remaining individuals who have phenotypes but no genotypes. We then select for sequencing the individual who maximizes the sum of the imputability scores of the remaining individuals. This procedure repeats until no phenotyped individuals who can be sequenced remain. The procedure then repeats, this time selecting among unphenotyped individuals who have available DNA, and prioritizing them based on whether have their genotypes can inform the genotypes of additional phenotyped (but with no available DNA) individuals. 

Intuitively, the final list of individuals to sequence is based on whether they have a known phenotype themselves and on their ability to inform the genotypes of additional phenotyped individuals in the pedigree.

## License
This project is licensed under BSD License 2.0

## Contact
Jimmy Liu (jliu@nygenome.org), New York Genome Center


