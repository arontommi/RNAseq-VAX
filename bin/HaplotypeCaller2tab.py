#!/usr/bin/env python3


import click
from vcffunctions.vcfheader import *
from vcffunctions.vcf2tab import *

@click.command()
@click.option(
    '--inputfile', '-i',
    help='give it strelka a strelka file!',
)
@click.option(
    '--outputfile', '-o',
    help='give it outputfile!',
)
@click.option(
    '--sample', '-s',
    help='give it a sample',
)
def main(inputfile, outputfile, sample ):
    """
    simple run wrapper around around extracting columns from haplotypecaller runs, that have been annotated and return
    a csv.

    """
    # Extract header info from VCF
    data = get_header_dict(inputfile)
    csq_header = read_csq(data['INFO']['CSQ'])

    # create pandas df
    data = read_vcf(inputfile)
    # only keep variants that passed
    passvcf = data[data['FILTER'] == 'PASS']
    # expand various columns
    data = extract_columns(passvcf.iloc[:, -1], ':', list('GT:AD:DP:GQ:PL'.split(':')))
    AD_data = extract_columns_varLength(data["AD"], ',', ['DP_REF','DP_ALT'])
    passvcf_info = col2dict(passvcf['INFO'])
    passvcf.reset_index(drop=True, inplace=True)

    # merge various coulumns to one large dataframe
    df_merged = pd.concat([data, AD_data, passvcf_info, passvcf], axis=1).reindex(passvcf.index)

    # create  a seperate line for each alt and for each CSQ annotaton
    newdf = splitdataframe(df_merged,'ALT',',')
    repopulated = splitdataframe(newdf, 'CSQ', ',')

    # expand the csq column and annotate it.
    csq_df = expandandextract(repopulated['CSQ'], csq_header)

    # merge the csq and the final column
    all_extracted = pd.concat([csq_df, repopulated], axis=1)
    all_extracted['Sample'] = sample


    # filter based on at least 10 allele reads and a ref depth 20
    all_extracted = all_extracted[all_extracted["DP_ALT"] >= 10]
    all_extracted = all_extracted[all_extracted['DP_REF'] >= 20]
    all_extracted.to_csv(outputfile)
if __name__ == '__main__':
    main()
