# SeqDef

SeqDef is an R package for calculating Sequencing Deficiency scores based on phylogenetic position and data availability.

## Installation

You can install the development version from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("Andresdbp/SeqDef")
```

## Usage

```r
library(SeqDef)

SeqDef(tree = tree, df = df, data.col = 2, invert = TRUE, scale = TRUE, lambda = 1)
```
