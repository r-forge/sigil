# Ideas and todo notes for SIGIL course and `corpora` package

## corpora Package

### Version 0.6-1

- bug fix (2024-11-24): most functions now convert frequency data to floating-point format in order to avoid overflow and other issues (`bit64::integer64` from `fread()` is particularly catastrophic as it masequerades as `double`)
- implement **bootstrapping** for keyness measures, which should be reasonably feasible (as only text-level frequency tables are needed and can easily be aggregated)
- release as version **0.6-1** or **0.6-2**

### Version 0.7

- implement function `am.score()` for computing best-practice association measures from vectorised frequency signatures $f, f_1, f_2, N$ 

- supports both pre-defined AM as well as user-defined equations, using an implementation similar to `wordspace` (but written in pure R and without batch processing for large data sets)

- optional: bootstrapping for AM scores based on per-text frequency counts (which need to  be provided separately for $f, f_1, f_2$ and $N$); `am.bootstrap()` will provide efficient implementation using `data.table`, but this should be an optional dependency (or at least only loaded on demand)

- probably better to provide bootstrapping code as a separate publication (possibly for a _Corpus Linguistics_ or _ICAME_ talk?) or a package vignette (only built if `data.table` is available) to keep `corpora` free of unnecessary dependencies

### Future plans

- integrate R6 class and support functions for Geometric Multivariate Analysis (GMA) into the `corpora` package
- first need some additional testing and usage to finalise the API

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
