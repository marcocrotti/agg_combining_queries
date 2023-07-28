process FIND_CHUNK {
    /*
    Find the genomic and annotated aggV2 vcf chunk
    */

    publishDir (
        "${params.outdir}",
        mode: 'copy',
        saveAs: { file -> 'find_chunk/' + file },
        pattern: "{*_files.txt}",
        enabled: params.publish_all
    )


    input:
    path(my_bed)
    path(aggv2_bed)

    output:
    path "geno_files.txt", emit: geno_vcf_list
    path "anno_files.txt", emit: vep_vcf_list
    path "versions.yml", emit : ch_versions_find_chunk

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
    echo "$gene,$gvcf,$gvcf_index" >> geno_files.txt
    echo "$gene,$avcf,$avcf_index" >> anno_files.txt
    done < !{my_bed}

    cat <<-EOF > versions.yml
    "${task.process}":
      bedtools: \$( bedtools --version | head -n1 | cut -d' ' -f2 )
    EOF

    '''

}