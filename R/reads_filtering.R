
# I used a one liner bash script from biostars: https://www.biostars.org/p/139006/#332214
# for i in `ls *.fastq.gz`; do echo "$i"; echo $($(zcat ${i} | wc -l)/4|bc); done
# It was modified from a question from biostars
# Read file and read
seqs <- read.delim("output/reads.txt", sep = " ",header = FALSE)
pdf("Figures/read_sequences2.pdf")
barplot(sort(seqs$V2), ylim = c(0, 10 ^5))
barplot(log10(sort(seqs$V2)), ylim = c(0, 5))
dev.off()

# Filter by the agreed threshold
s2 <- seqs[(seqs$V2 > 3000), ]
fs <- file.path("data/fastq_ASV/", s2$V1)

# Write a single file per line
write.table(s2[, 1, drop = FALSE], "output/samples_3000.csv", row.names = FALSE, col.names = FALSE,
            quote = FALSE)

# Scrip used to copy the files
while read p; do
  scp $(echo "../data/fastq_ASV/$p lrevilla@servidor-ciberehd.upc.es:/srv/juanjo/lrevilla/BCN_3000/")
  done < samples_3000.csv
