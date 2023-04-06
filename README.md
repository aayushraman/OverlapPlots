# README

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons License" style="border-width:0"/></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

## Why use Overlap plots?

Recent studies have suggested that genes longer than 100 kb are more likely to be misregulated in neurological diseases associated with synaptic dysfunction, such as autism and Rett syndrome (Gabel et al. 2015, Boxer et al. 2020). These length-dependent transcriptional changes are modest in MeCP2-mutant samples, but, given the low sensitivity of high-throughput transcriptome profiling technology, here we re-evaluate the statistical significance of these results. We find that the apparent length-dependent trends previously observed in MeCP2 microarray and RNA-sequencing datasets disappear after estimating baseline variability from randomized control samples. Overlap plots helps compare the KO/WT and WT/WT variations.

## What are the changes from the previous version of the OverlapPlot (Raman et al., *Nature Comm* 2018)?

We have included functions for generating two groups of WT clusters to assess WT/WT variation. These are: 

1. original *k*-means
2. variation of *k*-means function to generate clusters of equal sizes since the original will not always outputs clusters of equal sizes.

## Some examples of Overlap Plots for Boxer et al., *Molecular Cell* 2020 dataset

**1. Overlap plot for KO/WT whole cell dataset**

<img src='dat/ex_overlap_plots/KO-WT_whole-cell.png' width='300'>

**2. Overlap plot for KO/WT nuclear dataset**

<img src='dat/ex_overlap_plots/KO-WT_nuclear.png' width='300'>

**3. Overlap plot for KO/WT chromatin dataset**

<img src='dat/ex_overlap_plots/KO-WT_chromatin.png' width='300'>

## References

1.  **AT Raman**†, AE Pohodich†, YW Wan, HK Yalamanchili, HY Zoghbi, Z Liu. Apparent bias towards long gene misregulation in MeCP2 syndromes disappears after controlling for baseline variations. [*Nature Communications* (2018)](https://www.nature.com/articles/s41467-018-05627-1) (PMID: 30104565)

2.  **AT Raman**. A research parasite's perspective on establishing a baseline to avoid errors in secondary analyses. [*GigaScience* (2021)](https://academic.oup.com/gigascience/article/10/3/giab015/6168809) (PMID: 33710326)

## Inquiries

Please add issues if you have any questions/comments regarding the "overlap plots" or please contact directly at [aayushraman09\@gmail.com](mailto:aayushraman09@gmail.com).

## Contributors

1.  Ayush T Raman [original author]
2.  Daniel Palacios [currently working on the bioconductor package]
3.  Zhandong Liu [supervised the work]
