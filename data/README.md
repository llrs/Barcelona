The data of P4 folder (from plate 4) and other plates is coming from [Partek](http://partek.servidor-ciberehd.upc.es/), processed with Kraken.


They are processed using [Kraken-biom](https://github.com/smdabdoub/kraken-biom)
`kraken-biom *001-w000_S37_L001_R1_001*.fastq_output.txt --max D --fmt json` 
but it doesn't report a valid JSON as stated with :

```
if jq -e . >/dev/null 2>&1 <<<"$json_string"; then
    echo "Parsed JSON successfully and got something other than false/null"
else
    echo "Failed to parse JSON, or got false/null"
fi
```
to be then imported with pyloseq function import_biom


The Silva data from version 138 release used for imputation was [downloaded](https://doi.org/10.5281/zenodo.3986799 "Zendo site with DADA2 ready silva database")(silva_nr99_v138_train_set.fa.gz) on the 2020/12/03.
