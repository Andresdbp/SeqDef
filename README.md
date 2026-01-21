# SeqDef: Sequencing Deficiency & Prioritization

**SeqDef** is an R package designed to prioritize taxa for genomic sequencing initiatives (such as the Earth BioGenome Project). It quantifies the marginal genomic value of individual taxa by weighting the availability of sequence resources in related species using a distance-decay kernel.

Unlike static priority lists, SeqDef scores are dynamic: if a close relative of a target species is sequenced, the target's priority score decreases to reflect the shared information redundancy.

## Features

- **Phylogenetic Awareness:** Uses a distance-decay kernel to weigh information sharing across the tree.
- **Automated Calibration:** Includes two algorithms to automatically select the optimal decay parameter (`lambda`):
  - `auto_max`: Maximizes the variance of the statistic to ensure discriminatory power.
  - `by_genus`: Anchors the decay rate to the median evolutionary distance between congeneric species.
  - A manual value for `lambda` is also an option by providing a numeric value
- **Composite Prioritization:** Includes the function (`calc_priority`) to combine SeqDef scores with external variables (e.g., IUCN Red List status, economic value).

## Installation

You can install the development version of SeqDef from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Andresdbp/SeqDef")
```

## Quick Start

``` r
library(SeqDef)
library(ape)

# 1. Load your Tree and Data
# (Here we simulate a random tree for demonstration)
set.seed(42)
tree <- rtree(20)
# Data: 1 = Genome Available, 0 = Genome Missing
data <- data.frame(species = tree$tip.label, 
                   availability = sample(c(0,1), 20, replace = TRUE))

# 2. Calculate SeqDef Scores
# Uses "auto_max" to find the best lambda automatically
res <- SeqDef(tree, data, data.col = "availability", lambda = "auto_max")

# View top "Deficient" species (High SeqDef = High Priority)
head(sort(res$seqdef, decreasing = TRUE))

# 3. Calculate Composite Priority (e.g., SeqDef * 2^IUCN)
# Simulate IUCN scores (0 = LC, 4 = CR)
iucn_vals <- setNames(sample(0:4, 20, replace = TRUE), tree$tip.label)

final_priority <- calc_priority(res, iucn_vals, model = "exponential", base = 2)

# View final priority list
head(sort(final_priority, decreasing = TRUE))
```

## Methodology

SeqDef defines the Sequencing Deficiency of a focal taxon as the complement of its phylogenetically weighted data availability. The weight between any two species decays exponentially with evolutionary distance.The package solves the problem of "arbitrary parameter choice" by optimizing the decay parameter ($\lambda$) to fit the specific topology and data distribution of the user's clade.

## Citation

Citation
If you use SeqDef in your research, please cite:

[TBD]

