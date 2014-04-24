##
##  Prepare data sets for inclusion in corpora/SIGIL package
##

## BNC metadata
BNCmeta <- read.delim("tbl/bnc_metadata_utf8.tbl", quote="", fileEncoding="UTF-8", encoding="UTF-8")
save(BNCmeta, file="rda/BNCmeta.rda", compress="xz")

## text statistics for Brown and LOB corpora
BrownStats <- read.delim("tbl/brown.stats.txt", quote="")
save(BrownStats, file="rda/BrownStats.rda", compress="xz")

LOBStats <- read.delim("tbl/lob.stats.txt", quote="")
save(LOBStats, file="rda/LOBStats.rda", compress="xz")

## number of passives per genre in Brown and LOB corpora
BrownPassives <- read.csv("tbl/passives.brown.csv")
save(BrownPassives, file="rda/BrownPassives.rda", compress="xz")

LOBPassives <- read.csv("tbl/passives.lob.csv")
save(LOBPassives, file="rda/LOBPassives.rda", compress="xz")

## adjacent bigrams in the Brown corpus
BrownBigrams <- read.delim("tbl/brown_bigrams.tbl", quote="")
save(BrownBigrams, file="rda/BrownBigrams.rda", compress="xz")

## PP-verb collocations annotated by Brigitte Krenn
KrennPPV <- read.delim("tbl/krenn_pp_verb.tbl", quote="", fileEncoding="UTF-8", encoding="UTF-8")
save(KrennPPV, file="rda/KrennPPV.rda", compress="xz")

## passive counts for each text in the Brown and Lob corpora
BrownLOBPassives <- read.delim("tbl/brown_lob_passives.tbl", quote="")
save(BrownLOBPassives, file="rda/BrownLOBPassives.rda", compress="xz")

## -- tbl/bigrams.100k.tfl has to be loaded with zipfR package, so don't include in SIGIL


## make ZIP archive containing all data files
zip.name <- "sigil_datasets.zip"
if (file.exists(zip.name)) file.remove(zip.name)
system2("zip", c("-r", zip.name, "tbl/*"))
