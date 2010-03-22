##
##  Prepare data sets for inclusion in corpora package
##

# -- BNC metadata
BNCmeta <- read.delim("bnc_metadata_utf8.tbl", quote="", fileEncoding="UTF-8", encoding="UTF-8")
save(BNCmeta, file="bnc_metadata.rda")
