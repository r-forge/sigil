# Ideas and todo notes for SIGIL course and `corpora` package

## corpora Package

### Version 0.7

- new release for talks on boostrapping at ICAME 2025 and Corpus Linguistics 2025 $\to$ should be on CRAN by 15 June 2025 (so submit around 9 June)

- optional: implement **bootstrapping** for keyness measures, which should be reasonably feasible (as only text-level frequency tables are needed and can easily be aggregated)

- optional: implement bootstrapping for AM scores based on per-text frequency counts with function `am.bootstrap()` 
  - efficient implementation might need `data.table` or `tidyverse` for aggregation, but `corpora` shouldn't depend on these packages
  - if there is no good base R implementation, better not to include in package to keep `corpora` free of unnecessary dependencies
  
- perhaps better to provide bootstrapping code in online repo (OSF or on R-Forge homepage) accompanying the ICAME / Corpus Linguistics talks
  - can use any suitable packages for efficiency and convenience; might also show both `tidyverse` and `data.table` approach
  - could still be integrated as `corpora` functions in a future release

### Future plans

- none at the moment

## SIGIL Course

### General

- make HTML renderings of worked examples (`Rmd`) available on SIGIL homepage

### Unit #2

- integrate LaTeX slides from FAU course instead of old Keynote presentation
- check that all important content is included

### Unit #3

- should be revised based on recent statistics courses at FAU, including worked example
- show Anscombe's quartet as a warning on trusting summary statistic
  - see http://en.wikipedia.org/wiki/Anscombe's_quartet
  - included in base R as data set `anscombe` 

### Unit #4

- should be reworked based on recent (simpler) teaching materials on keywords and collocations
- possibly use existing PowerPoint slides as an interim solution

### Unit #6

- create this unit based on material from recent statistics courses at FAU

### Unit #7

- expand slides on multivariate data analysis
- ideally including overview of Biber-style MDA and the GMA approach; as an interim solution, the PowerPoint presentation can be used

### Unit #8

- redo illustration plots for passives example
  - R code already imported into SVN, revise as necessary
- re-structure and improve discussion of non-randomness
- develop exercise using the passives data (or perhaps noun frequencies)
  - LM / GLM regression
  - either include register/dsm predictors in SIGIL package
  - or use data set from BNC with more comprehensive metadata
- ideally to be recreated as a LaTex presentation, but may be too much unnecessary work
