#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """\

    C O M P L E X  V A R I A N T
    ============================
    pipeline version  : 1.0 alpha
    input_bed         : ${input_bed}
    agg_chunks_bed    : ${params.agg_chunks_bed}
    severity_scale    : ${params.severity_scale}
    severity          : ${params.severity}
    include_exclude   : ${include_exclude}
    expression        : ${expression}
    format            : ${format}
    output directory  : ${params.outdir}
    command           : ${workflow.commandLine}
    """.stripIndent()

/*---------------------
  Start pipeline
 ----------------------*/

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

include { FIND_CHUNK } from "../modules/local/find_chunk/find_chunk.nf"
include { EXTRACT_VARIANT_VEP_SEVERITY_SCALE } from "../modules/local/extract_variant_vep_severity_scale/extract_variant_vep_severity_scale.nf"
include { EXTRACT_VARIANT_VEP } from "../modules/local/extract_variant_vep/extract_variant_vep.nf"
include { INTERSECT_ANNOTATION_GENOTYPE_VCF } from "../modules/local/intersect_annotation_genotype_vcf/intersect_annotation_genotype_vcf.nf"
include { FIND_SAMPLES } from "../modules/local/find_samples/find_samples.nf"
include { SUMMARISE_OUTPUT } from "../modules/local/summarise_output/summarise_output.nf"

workflow {

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