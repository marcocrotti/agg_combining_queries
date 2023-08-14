#!/bin/bash

#BSUB -P Bio
#BSUB -q inter
#BSUB -J agg_combine_query
#BSUB -o agg_combine_query_%J.stdout
#BSUB -e agg_combine_query_%J.stderr

module purge
module load bio/nextflow/22.10.5 tools/singularity/3.8.3

# Location of the Small Variant workflow
agg_combine_query='/pgen_int_work/BRS/mcrotti/tinkering/agg_combining_queries'

nextflow run "${agg_combine_query}"/main.nf \
    --publish_all true \
    -profile cluster \
    -resume