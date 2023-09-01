#!/usr/bin/env python3
# Adapted from: https://github.com/nf-core/modules/tree/master/modules/nf-core/custom/dumpsoftwareversions

import click
import yaml


@click.command()
@click.option('--versions', help='YAML versions files.')
@click.option('--command', help='Command line call.')
@click.option('--out', help='Name of the output file.')
def dump_versions(versions, command, out):
    """Load all version files and generate merged output."""

    with open(versions) as f:
        versions_by_process = yaml.load(f, Loader=yaml.BaseLoader)

    versions_by_module = {}
    for process, process_versions in versions_by_process.items():
        module = process.split(':')[-1]
        try:
            if versions_by_module[module] != process_versions:
                raise AssertionError(
                    'We assume that software versions are the same between all modules. '
                    'If you see this error-message it means you discovered an edge-case '
                    'and should open an issue in nf-core/tools. '
                )
        except KeyError:
            versions_by_module[module] = process_versions

    with open(out, 'w') as f:
        f.write(command)
        yaml.dump(versions_by_module, f, default_flow_style=False)


if __name__ == '__main__':
    dump_versions()