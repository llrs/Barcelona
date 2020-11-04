# Read files to copy ####
meta <- readRDS("data_out/refined_meta.RDS")
samples <- meta$Name
stopifnot(length(samples) == 128)

# Copy from the disk to the computer ####
out <- lapply(samples, function(x) {
  loc <- list.files(path = "/media/lrevilla/Elements/SEQ/Michigan/Complete/",
                    recursive = TRUE, full.names = TRUE,
                    pattern = paste0("^", x, ".*.fastq.gz"))
  dest <- "/home/lrevilla/Documents/projects/design_ngs/data/fastq_ASV/"
  file.copy(from = loc, to = dest, copy.date = TRUE, overwrite = FALSE)
})
summary(unlist(out))
# Remove some files I souldn't copy ####
dest <- "/home/lrevilla/Documents/projects/design_ngs/data/fastq_ASV"
keep <- lapply(samples, function(x){
  list.files(path = dest, full.names = TRUE, pattern = x)})
keep <- unlist(keep)
all_files <- list.files(path = dest, full.names = TRUE)
remove_files <- all_files[!all_files %in% keep]
file.remove(remove_files)

remaining <- list.files(path = dest, full.names = FALSE)
stopifnot(length(remaining) == 2*length(samples))
