# Barcelona

This project is to provide the batch groups for the sequencing at Michigan
and integrate microbiome and transcriptome

The data is mangled in preprocessing.R, which later was used to design the batch.


# Options considered

Some of this options were from a question I [posted online](https://bioinformatics.stackexchange.com/q/4765/48)

- OSAT.R: Uses the [OSAT library](https://bioconductor.org/packages/OSAT)
- github_notebook.R: Uses a suggested notebook.
- own_algorithm.R: uses the  [experDesign](https://github.com/llrs/experDesign) package I specifically created for this

I went with my own_algorithm solution

# Reality check
 
After checking the samples it came out that the preprocessing was wrong. 
Aida send fewer samples (still for 4 batches).

TODO: Check the preprocessing (Not so important but would be nice)
TODO: Check the final batches and effect of the concentrations given
TODO: Prepare the names and relation with the databases
The data is in dataDelivered, with the exception of amplicon_names.xlsx which is a copy 
paste of AmmpliconProject...FinalProject.xlsx. (see [this](https://github.com/tidyverse/readxl/issues/513) problem)

The last two points are in post-sequencing.R
