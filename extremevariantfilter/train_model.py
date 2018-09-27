#!/home/eanderson/Virtual_Envs/Benchmark_Filtering/bin/python
"""
Usage:
    train_model.py (--true-pos STR) (--false-pos STR) (--type STR) [--out STR] [--njobs INT] [--verbose]

Description:
    Train a model to be saved and used with VCFs.

Arguments:
    --true-pos STR          Path to true-positive VCF from VCFeval or comma-seperated list of paths
    --false-pos STR         Path to false-positive VCF from VCFeval or comma-seperated list of paths
    --type STR              SNP or INDEL

Options:
    -o, --out <STR>                 Outfile name for writing model [default: (type).filter.pickle.dat]
    -n, --njobs <INT>               Number of threads to run in parallel [default: 2]
    -h, --help                      Show this help message and exit.
    -v, --version                   Show version and exit.
    --verbose                       Log output

Examples:
    python train_table.py --true-pos <path/to/tp/vcf> --false-pos <path/to/fp/vcf> --type [SNP, INDEL] --njobs 20
"""

# Load libraries
import pandas as pd
import numpy as np
import svf_leave_one_out_snp_par_one as svf
import confusion_matrix as cm
import pickle
import warnings
from xgboost import XGBClassifier
from sklearn.preprocessing import LabelBinarizer
from multiprocessing import Pool
from contextlib import closing
from docopt import docopt


def get_options():
    args = docopt(__doc__, version='1.0')
    verbose = args['--verbose']

    # Read training data
    tp_vcf = args['--true-pos']
    fp_vcf = args['--false-pos']
    poly = args['--type']
    njobs = int(args['--njobs'])
    outname = args['--out']
    if outname == "(type).filter.pickle.dat":
        outname = poly + '.filter.pickle.dat'
    cm.check_type(poly)

    return tp_vcf, fp_vcf, poly, njobs, outname


def Get_Name_Mod(variant_table, poly):
    filename = variant_table.split('/').pop()
    basename = '.'.join(filename.split('.')[0:-1])
    outname = basename + '.pickle.dat'
    return filename, outname


def Open_VCF(vcf_path):
    header=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'CALLS']
    #vcf = pd.read_csv(vcf_path, delimiter="\t", comment="#", names=header, compression="gzip")
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


def Make_Table(vcf, label):
    var_vcf = Open_VCF(vcf)
    info_fields = pd.DataFrame(list(var_vcf['INFO'].apply(Split_Info))).fillna(0.)
    info_fields = info_fields[['QD', 'MQ', 'FS', 'MQRankSum', 'ReadPosRankSum', 'SOR']]
    calls = Get_Calls_Info(var_vcf)
    info_fields = pd.concat([info_fields, calls], axis=1)
    info_fields['label'] = label
    return info_fields


def Get_Training_Table(tp_vcf, fp_vcf):
    tp_snp_vcf = Make_Table(tp_vcf, 1)
    fp_snp_vcf = Make_Table(fp_vcf, 0)
    full_vcf = tp_snp_vcf.append(fp_snp_vcf)
    array = full_vcf.values
    X = array[:,0:11]
    Y = array[:,11]
    return X, Y.astype(int)


def Get_Training_Tables(tp_fp_tup):
    tp_vcf, fp_vcf = tp_fp_tup
    tp_snp_vcf = Make_Table(tp_vcf, 1)
    fp_snp_vcf = Make_Table(fp_vcf, 0)
    full_vcf = tp_snp_vcf.append(fp_snp_vcf)
    array = full_vcf.values
    X = array[:,0:11]
    Y = array[:,11]
    return X, Y.astype(int)


def Check_VCF_Paths(tp_vcf_path, fp_vcf_path):
    if "," in tp_vcf_path or "," in fp_vcf_path:
        if len(tp_vcf_path.split(",")) != len(fp_vcf_path.split(",")):
            raise ValueError('Unequal number of True and False VCFs supplied')

    return zip(tp_vcf_path.split(","), fp_vcf_path.split(","))


if __name__ == "__main__":
    warnings.filterwarnings('ignore',category=DeprecationWarning)

    tp_vcf, fp_vcf, poly, njobs, outname = get_options()
    all_vcf = Check_VCF_Paths(tp_vcf, fp_vcf)
    with closing(Pool(processes=njobs)) as pool:
        results = pool.map(Get_Training_Tables, all_vcf)
        pool.terminate()

    X, Y = zip(*results)
    X_all = np.concatenate(X)
    Y_all = np.concatenate(Y)

    model = cm.Build_Model(poly, njobs)
    print("Training {} on {}").format(model[0], all_vcf)
    model[1].fit(X_all, Y_all)
    with open(outname, "wb") as out:
        pickle.dump(model[1], out)
