/*---------------------
Intersect functional annotation VCF with genotype VCF
----------------------*/

process INTERSECT_ANNOTATION_GENOTYPE_VCF {
    /*
    Intersect functional annotation VCF with genotype VCF
    */

    input:
    tuple val(gene), path(gvcf), path(gvcf_index), path(avcf_subset), path(avcf_subset_index)

    output:
    tuple val(gene), path("${gene}_intersect/0000.vcf.gz"), path("${gene}_intersect/0000.vcf.gz.tbi"), emit : intersect_vcf
    path "versions.yml", emit : ch_versions_intersect_annotation_genotype_vcf

    script:

    """
    bcftools isec -i ${params.expression} -e- -p ${gene}_intersect -n=2 -O z ${gvcf} ${avcf_subset}

    cat <<-EOF > versions.yml
    "${task.process}":
      bcftools: \$( bcftools --version | head -n1 | cut -d' ' -f2 )
    EOF
    """

}
