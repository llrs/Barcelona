tab <- read.delim("data/P4/Partek_Michigan3_Kraken_Classified_species.txt")
counts <- tab[, -1]
seqs <- colSums(counts)
not_empty <- function(x) {
  sum(x != 0)
}
sp <- apply(counts, 2, not_empty)
ctrls <- grep("X500_", colnames(counts))

sp_ctrls <- apply(counts[, ctrls], 1, not_empty)

seqs_b <- colSums(counts[!sp_ctrls, ])
sp_b <- apply(counts[!sp_ctrls, ], 2, not_empty)
p4 <- cbind(sp, seqs, sp_b, seqs_b)

# Points below these are samples with as much sequences and species as the blanks!!
plot(log10(sp), log10(seqs))
abline(v = log10(max(sp[ctrls])))
abline(h = log10(max(seqs[ctrls])))
p4[ctrls, ] # Blanks with 34 and 19 species; they overlap in 10 species
