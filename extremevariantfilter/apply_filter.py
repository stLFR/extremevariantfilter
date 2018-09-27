#!/home/eanderson/Virtual_Envs/Benchmark_Filtering/bin/python
"""
Usage:
    apply_filter.py (--vcf STR) (--snp-model STR) (--indel-model STR) [--verbose]

Description:
    Apply models from train_model.py to a vcf

Arguments:
    --vcf STR                     VCF to be filtered
    --snp-model STR               Model for applying to SNPs
    --indel-model INT             Model for applying to InDels

Options:
    -h, --help                      Show this help message and exit.
    -v, --version                   Show version and exit.
    --verbose                       Log output

Examples:
    python apply_filter.py --vcf <table> --snp-model <snp.pickle.dat> --indel-model <indel.pickle.dat>
"""

import pandas as pd
import numpy as np
import pickle
import warnings
from docopt import docopt


def get_options():
    args = docopt(__doc__, version='1.0')
    verbose = args['--verbose']

    # Read training data
    vcf = args['--vcf']
    snp_mod = args['--snp-model']
    ind_mod = args['--indel-model']

    return vcf, snp_mod, ind_mod


def Get_Header(vcf_path):
    with open(vcf_path, 'r') as vcf:
        header = []
        newline = vcf.readline()
        while newline.startswith('#'):
            if 'FILTER' in newline or "CHROM\tPOS" in newline and filter_written == False:
                header.append('##FILTER=<ID=XGB_SNP,Description="Likely FP SNP as determined by loaded model">\n')
                header.append('##FILTER=<ID=XGB_IND,Description="Likely FP InDel as determined by loaded model">\n')
                filter_written = True
            header.append(newline)
            newline = vcf.readline()
    return header


def Open_VCF(vcf_path):
    header=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'CALLS']
    vcf = pd.read_csv(vcf_path, delimiter="\t", comment="#", names=header)
    return vcf


def Split_Info(info):
    fields = ['QD=', 'MQ=', 'MQRankSum=', 'FS=', 'ReadPosRankSum=', 'SOR=']
    parts = dict(part.split('=') for part in info.split(';') if any(field in part for field in fields))
    return parts


def Get_Calls_Info(vcf):
    call_fields = vcf['CALLS'].str.split(":", expand=True)
    call_fields.columns = vcf['FORMAT'][0].split(":")
    GTS = pd.get_dummies(call_fields['GT'])
    AD = call_fields['AD'].str.split(',', expand=True)
    AD.columns = ['RefD', 'AltD', 'AltAltD']
    AD = AD.drop('AltAltD', axis=1)
    AD['RefD'] = pd.to_numeric(AD['RefD'])
    AD['AltD'] = pd.to_numeric(AD['AltD'])
    AD['RDper'] = (AD['RefD']/(AD['RefD'] + AD['AltD']))
    AD['ADrat'] = (AD['AltD']/(AD['RefD'] + .1))
    calls = pd.concat([GTS['0/1'], AD], axis=1)
    return calls


def Check_SNP(vcf):
    if "," in vcf['ALT']:
        if len(vcf['REF']) == 1 and (len(vcf['ALT'].split(',')[0]) == 1 or \
                                     len(vcf['ALT'].split(',')[1]) == 1):
            return 1
        else:
            return 0
    elif len(vcf['REF']) == 1 and len(vcf['ALT']) == 1:
        return 1
    else:
        return 0


def Predict_Var(vcf, snp_mdl, ind_mdl):
    params = vcf.iloc[0:11].values
    if vcf['Is_SNP'] == 1:
        return int(snp_mdl.predict(params[None, :]))
    else:
        return int(ind_mdl.predict(params[None, :]))


def Add_Filter(vcf):
    if vcf['Is_SNP'] == 1 and vcf['Predict'] == 0:
        return "XGB_SNP"
    elif vcf['Is_SNP'] == 0 and vcf['Predict'] == 0:
        return "XGB_IND"
    else:
        return "."


def Get_Name(path):
    filename = path.split('/').pop()
    basename = '.'.join(filename.split('.')[0:-1])
    outname = basename + '.filter.vcf'
    return outname


def Write_VCF(table, header, outname):
    with open(outname, 'w') as vcf_out:
        for line in header:
            vcf_out.write("%s" % line)
    table.to_csv(outname, sep = '\t', header = False, mode = 'a', index = False)


if __name__ == "__main__":
    warnings.filterwarnings('ignore',category=DeprecationWarning)

    vcf_path, snp_mod, ind_mod = get_options()
    header = Get_Header(vcf_path)
    vcf = Open_VCF(vcf_path)

    with open(snp_mod, "rb") as snp_m:
        snp_mdl = pickle.load(snp_m)
    with open(ind_mod, "rb") as ind_m:
        ind_mdl = pickle.load(ind_m)

    info_fields = pd.DataFrame(list(vcf['INFO'].apply(Split_Info))).fillna(0.)
    info_fields = info_fields[['QD', 'MQ', 'FS', 'MQRankSum', 'ReadPosRankSum', 'SOR']]
    calls = Get_Calls_Info(vcf)
    info_fields = pd.concat([info_fields, calls], axis=1)
    info_fields['Is_SNP'] = vcf.apply(Check_SNP, axis=1)
    info_fields['Predict'] = info_fields.apply(Predict_Var, axis=1, args=(snp_mdl, ind_mdl))
    vcf['FILTER'] = info_fields.apply(Add_Filter, axis=1)

    Write_VCF(vcf, header, Get_Name(vcf_path))
