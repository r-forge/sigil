##
##  Prepare data sets for inclusion in corpora/SIGIL package
##

## BNC metadata
BNCmeta <- read.delim("tbl/bnc_metadata_utf8.tbl", quote="", fileEncoding="UTF-8", encoding="UTF-8")
BNCmeta$title <- as.character(BNCmeta$title)
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
BrownBigrams <- transform(BrownBigrams, word1=as.character(word1), word2=as.character(word2))
save(BrownBigrams, file="rda/BrownBigrams.rda", compress="xz")

## PP-verb collocations annotated by Brigitte Krenn
KrennPPV <- read.delim("tbl/krenn_pp_verb.tbl", quote="", fileEncoding="UTF-8", encoding="UTF-8", stringsAsFactors=FALSE)
save(KrennPPV, file="rda/KrennPPV.rda", compress="xz")

## passive counts for each text in the Brown and Lob corpora
BrownLOBPassives <- read.delim("tbl/brown_lob_passives.tbl", quote="")
save(BrownLOBPassives, file="rda/BrownLOBPassives.rda", compress="xz")

## -- tbl/bigrams.100k.tfl has to be loaded with zipfR package, so don't include in SIGIL

## BNC sample data for frequency comparison and collocation analysis
BNCInChargeOf <- read.delim("tbl/bnc_in_charge_of.tbl", quote="", stringsAsFactors=FALSE)
save(BNCInChargeOf, file="rda/BNCInChargeOf.rda", compress="xz")

BNCcomparison <- read.delim("tbl/bnc_comparison.tbl", quote="", stringsAsFactors=FALSE)
save(BNCcomparison, file="rda/BNCcomparison.rda", compress="xz")

BNCdomains <- read.delim("tbl/bnc_domains.tbl", quote="")
save(BNCdomains, file="rda/BNCdomains.rda", compress="xz")



## make ZIP archive containing all data files
zip.name <- "sigil_datasets.zip"
if (file.exists(zip.name)) stopifnot(file.remove(zip.name))
system2("zip", c("-r", zip.name, "tbl/*"))
