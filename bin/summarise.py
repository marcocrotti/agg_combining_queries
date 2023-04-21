"""
Script for SUMMARISE_OUTPUT process
"""
import sys
import os.path
import pandas as pd

if os.path.isfile(sys.argv[1]):
    input_file = sys.argv[1]
    par_var_tbl = pd.read_csv(input_file, sep='\t', header=None)

    par_tbl = par_var_tbl.groupby([0]).size().reset_index(
        name='count').sort_values(
        by='count', ascending=False).rename(
        columns = {0:'platekey'})

    par_var_tbl['variant'] = par_var_tbl[1].map(str) + '_' + par_var_tbl[
        2].map(str) + '_' + par_var_tbl[3].map(str) + '_' + par_var_tbl[4].map(str)

    var_tbl = par_var_tbl.groupby(['variant']).size().reset_index(
        name='count').sort_values(by='count', ascending=False)

    out_file_name_platekey = os.path.splitext(os.path.basename(
        input_file))[0].replace("_results", "_platekey_summary") + ".tsv"

    out_file_name_variant = os.path.splitext(os.path.basename(
        input_file))[0].replace("_results", "_variant_summary") + ".tsv"

    par_tbl.to_csv(out_file_name_platekey, sep='\t', index=False)
    var_tbl.to_csv(out_file_name_variant, sep='\t', index=False)

else:

    print('Error: input file does not exist')
