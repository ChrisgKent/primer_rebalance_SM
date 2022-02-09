# Currently Under active development. 
The full readme will be available soon 

# primer_rebalance_SM
This rebalances multiplex primer schemes as described in "Improving the evenness of SARS-CoV-2 genome coverage by titration of primer concentration"

# Installation 
To install
```
git clone https://github.com/ChrisgKent/primer_rebalance_SM
```
Activate/install the conda environment

```
cd primer_rebalance
conda env create -f primer_rebalance.yaml

conda activate nmaskgen
```

# Inputs
```
--config scheme=<file_to_scheme.bed> mapped_bed_dir=<dir> output_dir=results/
```
```scheme``` - Should provide a path to a .bed file, containing the locations of the primers used in the scheme. In the same format as PrimalScheme generates (ref, start, end, name, pool, strand, sequence)

```mapped_bed_dir``` - Point to a directory containing only the mapped .bam files. Such a the <name>.sorted.bam files from ncov-artic-nf pipeline 

```output_dir``` - The output directory, containing the results. Defaults to results


