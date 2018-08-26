#!/usr/bin/env python3

import io
import os
import pandas as pd
import click


def read_inputfile(path):
    lines = []
    with open(path, 'r') as f:
        for l in f:
            # select the CSQ header, assumes "##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format:" - improve.
            if l.startswith('##INFO=<ID=CSQ'):
                csq_holder = l.split(" ")[6].split('|')
            elif not l.startswith('#'):
                lines.append(l)
    # create table and drop the final columns
    return pd.read_table(io.StringIO(str.join(os.linesep, lines)), usecols=[*range(0, 30)]), csq_holder


def splitdataframe(df, target_column, separator):
    ''' df = dataframe to split,
    target_column = the column containing the values to split
    separator = the symbol used to perform the split
    returns: a dataframe with each entry for the target column separated, with each element moved into a new row.
    The values in the other columns are duplicated across the newly divided rows.
    '''
    def splitListToRows(row, row_accumulator, target_column, separator):
        split_row = row[target_column].split(separator)
        for s in split_row:
            new_row = row.to_dict()
            new_row[target_column] = s
            row_accumulator.append(new_row)
    new_rows = []
    df.apply(splitListToRows, axis=1, args=(
        new_rows, target_column, separator))
    new_df = pd.DataFrame(new_rows)
    return new_df


def expandandextract(df, csq):
    expandedcsq = splitdataframe(df, 'CSQ', ',')
    csqextracted = expandedcsq["CSQ"].str.split('|', expand=True)
    csqextracted.columns = csq
    alldata = pd.concat([expandedcsq, csqextracted], axis=1)
    alldataselect = alldata[['SampleName', 'CHROM', 'POS', 'REF', 'ALT', 'SYMBOL', 'Consequence', 'BIOTYPE', 'VARIANT_CLASS', 'Feature_type', 'IMPACT', 'FLAGS',
                             'NORMALT', 'NORMREF', 'TUMALT', 'TUMREF', 'TUMVAF', 'TUMVARFRACTION', 'Existing_variation', 'EUR_AF', 'gnomAD_NFE_AF', 'MAX_AF', 'ID', 'FILTER']]
    alldatarenamed = alldataselect.rename(columns={'SampleName': 'SAMPLE', 'Consequence': 'CONSEQUENCE', 'Feature_type': 'FEATURE_TYPE', 'NORMALT': 'ALT_READS_NORMAL', 'NORMREF': 'REF_READS_NORMAL', 'TUMALT': 'ALT_READS_TUMOR', 'TUMREF': 'REF_READS_TUMOR', 'TUMVAF': 'VARIANT_ALLELE_RATIO'})
    return alldataselect


@click.command()
@click.option(
    '--inputfile', '-i',
    help='seqtool melt file with header and sample name in first column',
)
@click.option(
    '--outputfile', '-o',
    help='outputfile!',
)
def main(inputfile, outputfile):
    testdb, csq = read_inputfile(inputfile)
    testdb_split = expandandextract(testdb, csq)
    testdb_split.to_csv('{}'.format(outputfile), sep='\t', index=False)


if __name__ == '__main__':
    main()
