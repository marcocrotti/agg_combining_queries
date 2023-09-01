process DUMP_VERSIONS {
    /*
    Merge the version.yml files.
    */

    publishDir(
        "${params.outdir}",
        mode: 'copy'
    )

    input:
    path(versions)

    output:
    path("software_versions.yml")

    script:
    """
    set -eoux pipefail

    dump_versions.py \
    --versions ${versions} \
    --command "command: ${workflow.commandLine}\n\n" \
    --out software_versions.yml
    """
}

