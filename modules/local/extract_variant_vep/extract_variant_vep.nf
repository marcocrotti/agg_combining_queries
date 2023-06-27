process EXTRACT_VARIANT_VEP {
    /*
    Extract variants in the functional annotation VCF 
    */

    input:
    tuple val(gene), path(avcf), path(avcf_index)

    output:
    tuple val(gene), path("${gene}_annotation.vcf.gz"), path("${gene}_annotation.vcf.gz.csi"), emit: annotation_vcf
    path "versions.yml", emit : ch_versions_extract_variant_vep

    script:

    """
    bcftools +split-vep -i 'SYMBOL="'"${gene}"'"' -c SYMBOL ${avcf} -O z -o ${gene}_annotation.vcf.gz
    bcftools index ${gene}_annotation.vcf.gz

    cat <<-EOF > versions.yml
    "${task.process}":
      bcftools: \$( bcftools --version | head -n1 | cut -d' ' -f2 )
    EOF
    """

}