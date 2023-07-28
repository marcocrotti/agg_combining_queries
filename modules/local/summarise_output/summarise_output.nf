process SUMMARISE_OUTPUT {
    /*
    Create summary tables
    */

    publishDir (
        "${params.outdir}",
        mode: 'copy',
        saveAs: { file -> 'summarise_output/' + file },
        pattern: "{*_summary.tsv}",
        enabled: params.publish_all
    )

    input:
    path(query_result)

    output:
    path("*_summary.tsv")
    path "versions.yml", emit : ch_versions_summarise_output

    script:

    """
    python3 summarise.py ${query_result}

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$( python --version | head -n1 | cut -d' ' -f2 )
    EOF

    """

}