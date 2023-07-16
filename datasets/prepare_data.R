##
##  Prepare data sets for inclusion in corpora/SIGIL package
##

library(plyr)


## BNC metadata
BNCmeta <- read.delim("tbl/bnc_metadata_utf8.tbl", quote="", fileEncoding="UTF-8", encoding="UTF-8")
BNCmeta$id <- as.character(BNCmeta$id)
BNCmeta$title <- as.character(BNCmeta$title)
save(BNCmeta, file="rda/BNCmeta.rda", compress="xz")


## VSS corpus in tabular format
fh <- pipe("cqp -f decode_vss.cqp", open="r")
VSS <- read.delim(fh, quote="", fileEncoding="ascii", colClasses="character")
close(fh)

VSS <- transform(VSS, sentence=as.integer(sub("^[a-z]+", "", sentence))) # recode sentence IDs as integers
VSS <- transform(VSS, story=factor(story, levels=unique(story))) # convert story as factor with original ordering of stories
save(VSS, file="rda/VSS.rda", compress="xz")


## text statistics for Brown and LOB corpora
BrownStats <- read.delim("tbl/brown.stats.txt", quote="")
save(BrownStats, file="rda/BrownStats.rda", compress="xz")

LOBStats <- read.delim("tbl/lob.stats.txt", quote="")
save(LOBStats, file="rda/LOBStats.rda", compress="xz")


## number of passives per genre in Brown and LOB corpora
BrownPassives <- read.csv("tbl/passives.brown.csv", stringsAsFactors=FALSE)
save(BrownPassives, file="rda/BrownPassives.rda", compress="xz")

LOBPassives <- read.csv("tbl/passives.lob.csv", stringsAsFactors=FALSE)
save(LOBPassives, file="rda/LOBPassives.rda", compress="xz")


## adjacent bigrams in the Brown corpus
BrownBigrams <- read.delim("tbl/brown_bigrams.tbl", quote="")
BrownBigrams <- transform(BrownBigrams, word1=as.character(word1), word2=as.character(word2))
BrownBigrams <- transform(BrownBigrams, O11=as.numeric(O11), O12=as.numeric(O12), O21=as.numeric(O21), O22=as.numeric(O22))
save(BrownBigrams, file="rda/BrownBigrams.rda", compress="xz")


## PP-verb collocations annotated by Brigitte Krenn
KrennPPV <- read.delim("tbl/krenn_pp_verb.tbl", quote="", fileEncoding="UTF-8", encoding="UTF-8", stringsAsFactors=FALSE)
save(KrennPPV, file="rda/KrennPPV.rda", compress="xz")


## lookup vector for genre labels (indexed by section code)
brown.genres <- structure(as.character(BrownPassives$name), names=as.character(BrownPassives$cat))

## passive counts for each text in the Brown and Lob corpora
BrownLOBPassives <- read.delim("tbl/brown_lob_passives.tbl", quote="", stringsAsFactors=FALSE)
BrownLOBPassives <- transform(BrownLOBPassives,
  genre = factor(genre, levels=brown.genres), # ensure genres are listed in "natural" order
  cat = factor(cat, levels=names(brown.genres)),
  lang = factor(lang, levels=c("AmE", "BrE"))
)
BrownLOBPassives <- droplevels(BrownLOBPassives) # remove genres not included in the table
save(BrownLOBPassives, file="rda/BrownLOBPassives.rda", compress="xz")


## passive counts wrt. VP tokens for each text in the Brown Family (courtesy of Gerold Schneider)
PassBFtokens <- read.delim("tbl/brownfam_vp_tokens.tbl.gz", stringsAsFactors=FALSE)
PassiveBrownFam <- ddply(PassBFtokens, .(id, corpus, section, genre, period, lang, n.words), function (X) table(X$voice))
PassiveBrownFam <- transform(PassiveBrownFam,
  corpus=factor(corpus, levels=c("BLOB", "Brown", "LOB", "Frown", "FLOB")),
  section=factor(section, levels=names(brown.genres)),
  genre=factor(genre, levels=brown.genres),
  period=factor(period),
  lang=factor(lang, levels=c("AmE", "BrE")),
  verbs=act+pass,
  p.pass=100*pass/(act+pass))
rownames(PassiveBrownFam) <- PassiveBrownFam$id
save(PassiveBrownFam, file="rda/PassiveBrownFam.rda", compress="xz")


## corresponding distributional features for advanced GLM analysis
load("raw/distfeat_brownfam.rda", verbose=TRUE)
text.id <- grep("^(brown|lob|frown|flob|blob)", rownames(BrownLReg10), value=TRUE, perl=TRUE) # text IDs we want to use
stopifnot(all(text.id %in% rownames(BrownLRegV10)))
stopifnot(all(text.id %in% rownames(BrownLTop10)))
DistFeatBrownFam <- data.frame(id=text.id, row.names=text.id, stringsAsFactors=FALSE)
DistFeatBrownFam <- cbind(DistFeatBrownFam, BrownLTop10[text.id, 1:9], BrownLReg10[text.id, 1:9], BrownLRegV10[text.id, 1:4])
tmp <- DistFeatBrownFam
for (i in 2:ncol(tmp)) tmp[, i] <- round(tmp[, i], 6)
write.table(tmp, file=gzfile("tbl/brownfam_distfeat.tbl.gz"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
save(DistFeatBrownFam, file="rda/DistFeatBrownFam.rda", compress="xz")

## -- tbl/bigrams.100k.tfl has to be loaded with zipfR package, so don't include in SIGIL


## BNC sample data for frequency comparison and collocation analysis
BNCInChargeOf <- read.delim("tbl/bnc_in_charge_of.tbl", quote="", stringsAsFactors=FALSE)
save(BNCInChargeOf, file="rda/BNCInChargeOf.rda", compress="xz")

BNCcomparison <- read.delim("tbl/bnc_comparison.tbl", quote="", stringsAsFactors=FALSE)
save(BNCcomparison, file="rda/BNCcomparison.rda", compress="xz")

BNCdomains <- read.delim("tbl/bnc_domains.tbl", quote="")
save(BNCdomains, file="rda/BNCdomains.rda", compress="xz")


## Biber features for texts in British National Corpus (from Gasthaus 2007)
BNCbiber <- read.delim("tbl/bnc_biber.tbl", quote="", check.names=FALSE)
ids <- BNCbiber[, 1]
BNCbiber <- as.matrix(BNCbiber[, -1])
rownames(BNCbiber) <- ids
colnames(BNCbiber) <- paste0("f", colnames(BNCbiber))
save(BNCbiber, file="rda/BNCbiber.rda", compress="xz")


## per-text frequency counts for some BNCweb queries (see prepare_data.tscript)
BNCqueries1 <- read.delim("raw/bnc_freqs_sentence.tbl", stringsAsFactors=FALSE)
BNCqueries2 <- read.delim("raw/bnc_freqs_word.tbl", stringsAsFactors=FALSE)
stopifnot(all(BNCqueries1$id == BNCqueries2$id))
BNCqueries <- cbind(BNCqueries1, subset(BNCqueries2, select=-id))
write.table(BNCqueries, file="tbl/bnc_queries.tbl", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
tmp <- read.delim("tbl/bnc_queries.tbl", stringsAsFactors=FALSE)
stopifnot(identical(tmp, BNCqueries)) # validate generated disk file
save(BNCqueries, file="rda/BNCqueries.rda", compress="xz")


## make ZIP archive containing all data files
if (FALSE) {
  ## should rather pack selected files manually
  zip.name <- "sigil_datasets.zip"
  if (file.exists(zip.name)) stopifnot(file.remove(zip.name))
  setwd("tbl")
  system2("zip", c("-r", paste0("../", zip.name), "*"))
  setwd("..")
}
