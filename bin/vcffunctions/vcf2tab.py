import io
import os
import pandas as pd
import click
from collections import defaultdict

def read_vcf(path):
    """

    :param path: path to vcf file
    :return: pandas dataframe with variants
    """
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_table(
        io.StringIO(str.join(os.linesep, lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str}
    ).rename(columns={'#CHROM': 'CHROM'})


def extract_columns(df_column, spliter, col_names):
    extracted_column = pd.DataFrame(df_column.str.split('{}'.format(spliter),
        n=len(df_column.head().str.get_dummies(sep='{}'.format(spliter)).columns), expand=True))
    extracted_column.columns = col_names
    extracted_column.reset_index(drop=True, inplace=True)
    return extracted_column


def col2dict(df_column):
    dictlist = []
    for l in df_column.str.split(';'):
        dictlist.append(dict((s.split('=') + [1])[:2] for s in l))
    df = pd.DataFrame(dictlist)
    df.reset_index(drop=True, inplace=True)
    return df

def expandandextract(df_csq_col, csq_holder):
    csqextracted = df_csq_col.str.split('|',expand=True)
    csqheader = csq_holder
    csqextracted.columns = csqheader
    return csqextracted




def splitdataframe(df,target_column,separator):
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
    df.apply(splitListToRows,axis=1,args = (new_rows,target_column,separator))
    new_df = pd.DataFrame(new_rows)
    return new_df

