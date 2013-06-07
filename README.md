cbcbSEQ: RNAseq analysis for UMD CBCB collaborators
====================================================
  
The purpose of this pipeline is to streamline the process for analyzing RNA-seq data
with potential batch effects. The pipeline includes 1) quantile normalization 2) log-transformation of counts 3) ComBat (location) batch correction 4) voom calculation of weights.

The functions in this package can be grouped into two main categories:

1. The functions used for assessing batch effects.
    a. `makeSVD`
    b. `pcRes`
    c.  `plotPC`
2. The functions for removing batch effect and computing weights for limma.
    a. `qNorm`
    b. `log2CPM`
    c.  `voomMod`
    d.  `combatMod`
    e.  `batchSEQ`
    
`batchSEQ` is the pipeline function. It combines `qNorm`, `log2CPM`, `voomMod`, and `combatMod`
into one step.

## Installation

```{r}
require(devtools)
install_github("cbcbSEQ", user="kokrah")
```