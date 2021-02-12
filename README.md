# Barcelona

This project is to provide the batch groups for the sequencing at Michigan
and integrate microbiome and transcriptome.

The data is mangled in preprocessing.R, which later was used to design the batch.


# Batch creation: options considered

Some of this options were from a question I [posted online](https://bioinformatics.stackexchange.com/q/4765/48)

- OSAT.R: Uses the [OSAT library](https://bioconductor.org/packages/OSAT)
- github_notebook.R: Uses a suggested notebook.
- own_algorithm.R: uses the  [experDesign](https://github.com/llrs/experDesign) package I specifically created for this

I went with my own_algorithm solution

# Done

The data is in dataDelivered, with the exception of amplicon_names.xlsx which is a copy 
paste of AmmpliconProject...FinalProject.xlsx. (see [this](https://github.com/tidyverse/readxl/issues/513) problem)

The last two points are in post-sequencing.R

Integration is done using RGCCA and our published method inteRmodel.
Steps:
1. pre-processing.R # Prepares data
2. modelling.R | models2.R | models3.R # Look for models, can be performed in parallel but models3 require some manual steps
3. models_report.R # Visualize together the results. 
4. botstrap_same_index.R # To bootstrap the samples and see if the models are accurate.
