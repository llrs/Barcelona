
# I used a one liner bash script to create the read.txt file something along $(zcat | wc -l) | bc
# It was modified from a question from biostars
# Read file and read
seqs <- read.delim("output/reads.txt", sep = " ",header = FALSE)
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
