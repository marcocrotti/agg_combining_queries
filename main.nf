#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*---------------------
  Start pipeline
 ----------------------*/

/*---------------------
  Find the genomic and annotated aggV2 vcf chunk 
 ----------------------*/

process FIND_CHUNK {
    
    publishDir "${params.outdir}/find_chunk_output", mode: 'copy'

    input:
    path(my_bed)
    path(aggv2_bed)

    output:
    path(geno_files), emit: geno_vcf_list
    path(anno_files), emit: vep_vcf_list

    shell:

    '''
    set -eoux pipefail
    
    while read -r line; do
    gene="$(echo "${line}"| awk '{print $4}')";
    printf "${line}" > ${gene}.bed;
    gvcf="$(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 10)";
    gvcf_index="$(echo $(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 10).csi)";
    avcf="$(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 11)";
    avcf_index="$(echo $(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 11).csi)";
    echo "$gene,$gvcf,$gvcf_index" >> geno_files
    echo "$gene,$avcf,$avcf_index" >> anno_files
    done < !{my_bed}

    '''

}

/*---------------------
  Extract variants in the functional annotation VCF 
 ----------------------*/
process EXTRACT_VARIANT_VEP_SEVERITY_SCALE {

    input:
    tuple val(gene), path(avcf), path(avcf_index)
    each path(severity_scale)

    output:
    tuple val(gene), path("${gene}_annotation.vcf.gz"), path("${gene}_annotation.vcf.gz.csi"), emit: annotation_vcf

    script:
    """
    bcftools +split-vep -i 'SYMBOL="'"${gene}"'"' -c SYMBOL -s worst:${params.severity}+ -S ${severity_scale} ${avcf} -O z -o ${gene}_annotation.vcf.gz
    bcftools index ${gene}_annotation.vcf.gz
    """

    }

process EXTRACT_VARIANT_VEP {

    input:
    tuple val(gene), path(avcf), path(avcf_index)

    output:
    tuple val(gene), path("${gene}_annotation.vcf.gz"), path("${gene}_annotation.vcf.gz.csi"), emit: annotation_vcf

    script:
    """
    bcftools +split-vep -i 'SYMBOL="'"${gene}"'"' -c SYMBOL ${avcf} -O z -o ${gene}_annotation.vcf.gz
    bcftools index ${gene}_annotation.vcf.gz
    """
    
}

/*---------------------
Intersect functional annotation VCF with genotype VCF
----------------------*/

process INTERSECT_ANNOTATION_GENOTYPE_VCF {

    input:
    tuple val(gene), path(gvcf), path(gvcf_index), path(avcf_subset), path(avcf_subset_index)

    output:
    tuple val(gene), path("${gene}_intersect/0000.vcf.gz"), path("${gene}_intersect/0000.vcf.gz.tbi")

    script:

    """
    bcftools isec -i ${params.expression} -e- -p ${gene}_intersect -n=2 -O z ${gvcf} ${avcf_subset}
    """

}

/*---------------------
Process results
----------------------*/

process FIND_SAMPLES {

    publishDir "${params.outdir}/final_output", mode: 'copy'

    input:
    tuple val(gene), path(int_vcf), path(int_vcf_index)

    output:
    path("${gene}_results.tsv")

    script:

    """
    bcftools query ${params.include_exclude} ${params.expression} -f ${params.format} ${int_vcf} > ${gene}_results.tsv
    """
}

/*----------------------
Create summary tables 
-----------------------*/

process SUMMARISE_OUTPUT {

    publishDir "${params.outdir}/final_output", mode: 'copy'

    input:
    path(query_result)

    output:
    path("*_summary.tsv") 

    script:

    """
    python3 summarise.py ${query_result}
    """

}

workflow {

    // user input bed file
    my_bed_ch = Channel
            .fromPath(params.input_bed)
            .ifEmpty { exit 1, "Cannot find input file : ${params.input_bed}" }

    // aggV2 bed chunks
    aggv2_bed_ch = Channel
            .fromPath(params.agg_chunks_bed)
            .ifEmpty { exit 1, "Cannot find input file : ${params.aggv2_chunks_bed}" }

    // VEP severity scale
    if (params.severity_scale) {
    severity_scale_ch = Channel
            .fromPath(params.severity_scale)
    }

    // Processes
    FIND_CHUNK(my_bed_ch, aggv2_bed_ch)	

    vep_vcf_ch = FIND_CHUNK.out.vep_vcf_list
        .splitCsv()
        .map {row -> [row[0], file(row[1]), file(row[2])] }

    geno_vcf_ch = FIND_CHUNK.out.geno_vcf_list
        .splitCsv()
        .map {row -> [row[0], file(row[1]), file(row[2])] }
    
    if (params.severity_scale != false) {
        EXTRACT_VARIANT_VEP_SEVERITY_SCALE(vep_vcf_ch, severity_scale_ch)
        intersect_input_ch = geno_vcf_ch
                                .join(EXTRACT_VARIANT_VEP_SEVERITY_SCALE.out.annotation_vcf)
    } else {
        EXTRACT_VARIANT_VEP(vep_vcf_ch)
        intersect_input_ch = geno_vcf_ch
                                .join(EXTRACT_VARIANT_VEP.out.annotation_vcf)
    }

	INTERSECT_ANNOTATION_GENOTYPE_VCF(intersect_input_ch)
	FIND_SAMPLES(INTERSECT_ANNOTATION_GENOTYPE_VCF.out)
	SUMMARISE_OUTPUT(FIND_SAMPLES.out.filter{ it.size()>0 })

}