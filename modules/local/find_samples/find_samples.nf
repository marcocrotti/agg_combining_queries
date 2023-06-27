process FIND_SAMPLES {
    /* 
    Process results
    */

    publishDir (
        "${params.outdir}",
        mode: 'copy',
        saveAs: { file -> 'find_samples/' + file },
        pattern: "{*_results.tsv}",
        enabled: params.publish_all
    )

    input:
    tuple val(gene), path(int_vcf), path(int_vcf_index)

    output:
    path("${gene}_results.tsv")

    script:

    """
    bcftools query ${params.include_exclude} ${params.expression} -f ${params.format} ${int_vcf} > ${gene}_results.tsv
    """
}