#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

def helpMessage() {
    println(
        """
        DESCRIPTION

        This workflow allows you to extract variants and samples that comply to both a set of genotype and functional annotation filters, by intersecting the genotype VCFs with the functional annotation VCFs.

        USAGE

        bsub < submit.sh

        PARAMETERS

        This is a subset of the common inputs a user may want to edit in the submission script.

        Query region parameters:

            --input_bed : This is a region file of your genes of interest. This must be a three or column tab-delimited file of chromosome, start, and stop (with an option fourth column of an identifier - i.e. a gene name). The file should have the .bed extension. 
            --expression : This parameter defines the bcftools filter of your query. See bcftools EXPRESSIONS for accepted filters.
            --format : This parameter defines the format of the query, see https://samtools.github.io/bcftools/bcftools.html#query for details. For the process to run, you should add the following fields '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\n]', but you can also specify additional fields after the initial list.

        List of aggregate VCFs:

            --agg_chunks_bed : This is the list of chunk names and full file paths to both the genotype and functional annotation VCFs for either aggV2 or aggCOVID.

        Optional parameters:

            --severity_scale : This file lists the severity of variants. Provide this file if interested only in variant with a specific consequence.
            --severity : With this parameter we choose the severity of variants we are interested in for our query. For example, if you want look only at missense variants or worse, the input value would be missense. Only use if the parameter severity_scale is set.


        """
    )
}

params.help = false
if (params.help){
    helpMessage()
    System.exit(0)
}

include { AGG_COMBINE_QUERIES } from "./workflows/agg_combine_queries.nf"

workflow {
    AGG_COMBINE_QUERIES()
}
