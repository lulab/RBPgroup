# RBPgroup

## Introduction

Two programs take a matrix containing non-negative values as an input, each column repersents
one CLIP-seq data, each row repersents one binding sites.
The optimal rank can be selected according to the results of `NMF.estimate.R`
The basis matrix and coefficient matrix as the main results can be used for further predicting 
RNA binding protein groups and potential functions for binding sites.

## Prerequisites
They are implemented in R and run under a UNIX/LINUX system.
To use the program, several R package should be installed appropriately:

* R (>3.0.0)
* R packages: argparse, NMF, bigmemory
* Python 2.7 is required for the `argparse` package to work.

We have tested our package using NMF version [0.17](https://cran.r-project.org/src/contrib/Archive/NMF/NMF_0.17.tar.gz) 
and [0.20.6](https://cran.r-project.org/src/contrib/NMF_0.20.6.tar.gz).

## Download
```bash
git clone https://github.com/lulab/RBPgroup.git
```

## Usage 
There are two core program named `NMF.estimate.R` and `NMF.main.R`.
To run them, it's like to run any other R program. If you don't specify any options, it will give the help options.

### Prepare the input matrix
The input file for NMF is a non-negative matrix, with features (RBPs) as columns and samples (RBP binding sites) as rows.
The input matrix file is a tab-delimited text file with column names specified in the first line and column names in the first
column. 
Usually, the `write.table()` function in R or `to_csv()` function in the pandas package in Python can be used to generate the input matrix file.
Each element of the matrix is the occupancy/binding affinity of an RBP on a RBP binding site, and is usually calculated as 
number of CLIP-seq reads in a RBP binding region normalized by number of reads in the same region from an RNA-seq experiment.

An data matrix used in our paper can be found in `data/example.zip`.
It is a zip archive and need to be extracted first.

The first few rows and columns of the matrix may look like this:
```
AGO1	AGO2	AGO3	AGO4	HNRNPA1	HNRNPA2B1
RE1000	0.14	0.04	0.16	0.12	0.01	0.21
RE10000	0.03	0.01	0.02	0.05	0.00	0.01
RE10001	0.05	0.05	0.03	0.02	0.00	0.18
RE10004	0.01	0.00	0.01	0.00	0.00    0.00
RE10009	0.02	0.00	0.00	0.04	0.00	0.03
RE1008	0.05	0.01	0.02	0.02	0.00	0.06
```

### Estimate the optimal rank
The script `NMF.estimate.R` is used to estimate the optimal rank for matrix factorization, using three quality measures:
cophenetic coefficient, dispersion coefficient and residuals. It takes a non-negative matrix file as input and outputs a
report of the metrics for rank selection.

Four NMF algorithms are included in the RBPgroup package: KL divergence (KL), Euclidean distance (euclidean),
KL divergence with orthogonality regularization (KL_ortho) and Euclidean distance with orthogonality regularization (euclidean_ortho).

The KL divergence is used by default because it leads to RBP clusters that are better supported by known PPI-networks
based on our experiments. The Euclidean distance is a common cost function for NMF, but tends to find a single cluster
with a large number of RBPs.

[iONMF](https://academic.oup.com/bioinformatics/article/32/10/1527/1742711/Orthogonal-matrix-factorization-enables) includes
an orthogonality regularization term that is useful for finding non-overlapping clusters. The weight on orthogonality
can be controlled by the hyper-parameter *alpha*. It is recommended to choose a value no larger than 10 for alpha.

```bash
bin/NMF.estimate.R -h
```
A help message will be shown:
```
usage: bin/NMF.estimate.R [-h] -i INPUT -o OUTPUT [-s RANK] [-e RANK]
                          [-m STRING] [-a NUMBER] [--seed NUMBER] [-n NUMBER]
                          [-p NUMBER]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input matrix
  -o OUTPUT, --output OUTPUT
                        output file prefix
  -s RANK, --start RANK
                        Number of start rank [default 1]
  -e RANK, --end RANK   Number of end rank
  -m STRING, --method STRING
                        the NMF algorithm to use. Should be one of
                        KL,euclidean,KL_ortho,euclidean_ortho. [default = KL]
  -a NUMBER, --alpha NUMBER
                        regularization factor for orthogonality of the
                        coefficient matrix [default = 10]
  --seed NUMBER         Seed for the random number generator
  -n NUMBER, --runs NUMBER
                        Number of runs to perform [default = 30]
  -p NUMBER, --processors NUMBER
                        Number of processors to use. This option is useful on
                        multicore *nix or Mac machine only, when performing
                        multiple runs (nrun > 1) [default 1]
```
I takes hours of time to perform NMF on a large matrix with many runs and ranks. It is recommended to 
specify the number of CPU cores through the `-p` option. Because NMF initializes the coefficient matrix and 
basis matrix with random values, specifying a number to the `--seed` option can make the NMF results reproducible.

output files:

* `*.txt`: this file contains values of three quality measures at different ranks.
* `*.pdf`: this file contains the line plots for three quality measures, heatmaps of consensus matrix for each ranks.

Run `NMF.estimate.R` on the example
```bash
bin/NMF.estimate.R -i data/example.mx -s 2 -e 5 -o output/estiRank/example.estiRank
```
Sample output files of the example can be found in `data/estiRank/`. The output files may be different each time
because random seed is selected based on the current time.

The most important quality measures in the output pdf file are cophenetic correlation coefficient (CPCC) and
dispersion coefficient. For both measures, a higher score indicates better results. The residuals are not very useful
for rank selection because they generally decrease as the rank increases. The ranks corresponding to local maximums 
of the cophenetic correlation coefficient and dispersion coefficient can be selected.

### Run NMF with the selected rank
The script `NMF.main.R` is used for non-negative matrix factorization with the selected rank. 
The input file is the same as the input file for `NMF.estimate.R`.
```bash
bin/NMF.main.R -h
```
A help message will be shown:
```
usage: bin/NMF.main.R [-h] -i INPUT -o OUTPUT -r RANK [--seed NUMBER]
                      [-n NUMBER] [-m STRING] [-a NUMBER] [-p NUMBER]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input matrix
  -o OUTPUT, --output OUTPUT
                        output file prefix
  -r RANK, --rank RANK  Number of rank
  --seed NUMBER         Seed for the random number generator
  -n NUMBER, --runs NUMBER
                        Number of runs to perform [default = 30]
  -m STRING, --method STRING
                        the NMF algorithm to use. Should be one of
                        KL,euclidean,KL_ortho,euclidean_ortho. [default = KL]
  -a NUMBER, --alpha NUMBER
                        regularization factor for orthogonality of the
                        coefficient matrix [default = 10]
  -p NUMBER, --processors NUMBER
                        Number of processors to use. This option is useful on
                        multicore *nix or Mac machine only, when performing
                        multiple runs (nrun > 1) [default = 1]
```

output files:

* `*.Rdata`: the raw data from non-negative matrix factorization 
* `*.pdf`: the heatmaps of basis matrix, coefficient matrix and consensus matrix
* `*.basis`: the basis matrix
* `*.coef`: the coefficient matrix
* `*.consensus`: the consensus matrix

Run `NMF.main.R` with rank 18 on the example:
```bash
bin/NMF.main.R -i data/example.txt -r 3 -n 100 -o output/main/example.3
```
The output files can be found in `data/main/example.18.*`.

### Extract the cluster components
The coefficient matrix contains the weights of each feature (RBP) in the clusters.
The script `NMF.assign_clusters.R` assign each feature to a cluster according to the weight of each feature
relative to the total weights of the features. It first normalizes each column (feature) of coefficient matrix
to sum to 1 by the sum of the column and assign a feature to a cluster if the normalized value is above a
predefined threshold. 

```
usage: bin/NMF.assign_clusters.R [-h] -i INPUT -o OUTPUT [-t THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        coefficient matrix
  -o OUTPUT, --output OUTPUT
                        output file for cluster components
  -t THRESHOLD, --threshold THRESHOLD
                        the threshold for coefficient matrix values to define
                        cluster components. [default = 0.2]
```
Run `NMF.assign_clusters.R` on the coefficient matrix:
```bash
bin/NMF.assign_clusters.R -i output/main/example.3.coef -o output/assign_clusters/example.3.assign_cluster.txt
```
A sample output file can be found in `data/assign_clusters/example.18.assign_cluster.txt`.


