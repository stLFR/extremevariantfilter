"""A machine learning approach to variant filtering using XGBoost

ExtremeVariantFilter is a machine learning based approach to variant filtering.
Provided models were trained on stLFR libraries from Complete Genomics.
These were generated from cell lines for HG001, HG002, HG003, HG004 and HG005.
This package provides functionality for applying existing models to VCFs,
as well as developing models based on personal or existing data.

Commands
--------
apply_filter (--vcf STR) (--snp-model STR) (--indel-model STR) [--verbose]
    Apply models generated by train_model to a VCF.

train_model (--true-pos STR) (--false-pos STR) (--type STR) [--out STR] [--njobs INT] [--verbose]
    Train a model to be saved and used to filter VCFs.

Functions
---------
Open_VCF(vcf_path)
    Given a vcf path returns a vcf loaded into memory as a DataFrame

Split_Info(info)
    Given a vcf info column, returns necessary info fields as a dict

Get_Calls_Info(vcf)
    Given a vcf returns encoded genotypes, reference and alternate
    allele depths, reference allele percentage and alternate
    allele ratio as a DataFrame

Check_Type(poly)
    Ensures input for polymorphism type is an accepted value

Check_VCF(vcf_path)
    Confirms that input VCF has only one calls field and
    a FORMAT field

Is_Gzipped(vcf_path)
    Checks to see whether a file has a .gz extension

Get_Header(vcf_path)
    Splits header from input vcf and stores it as a list

Check_SNP(vcf)
    Checks whether input variant is a SNP and returns a boolean

Predict_Var(vcf, snp_mdl, ind_mdl)
    Given a vcf, snp and indel models, returns a prediction as to
    whether or not the variant is real or now

Add_Filter(vcf)
    Updates the filter column value based on a variants predicted
    truth value

Get_Name(path)
    Given the path to the vcf, returns the root file name
    with '.filter.vcf' appended

Write_VCF(table, header, outname)
    Outputs the updated vcf with filter values added

Make_Table(vcf, label)
    Given a path for a true positive or false positive vcf and a label,
    loads the vcf into memory as a DataFrame and appends the label to
    the last column

Get_Training_Table(tp_vcf, fp_vcf)
    Given a true and false positive vcf paths, calls
    Make_Table() for both and then appends the
    DataFrames together

Get_Training_Tables(tp_fp_tup)
    Given a set of comma seperated true and false vcf paths,
    calls Make_Table for all paths and then appends the
    DataFrames together

Build_Model(poly, njobs)
    Returns an XGBoost classifier object with particular
    characteristics depending on whether it is being built
    to classify SNPs or InDels
"""

import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

import pandas as pd
import numpy as np
import gzip
import re
from xgboost import XGBClassifier
from sklearn.preprocessing import LabelBinarizer


# Common Functions


def Open_VCF(vcf_path):
    """Reads a vcf into memory

    Parameters
    ----------
    vcf_path : str
        File location of the vcf

    Returns
    -------
    DataFrame
        A DataFrame representation of the original vcf
    """

    header=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'CALLS']
    vcf = pd.read_csv(vcf_path, delimiter="\t", comment="#", names=header)
    return vcf


def Split_Info(info):
    """Splits necessary information out from the info vcf column

    Parameters
    ----------
    info : Series
        Info column from a vcf as a series

    Returns
    -------
    dict
        A dict of necessary fields and their values from the vcf
    """

    fields = ['QD=', 'MQ=', 'MQRankSum=', 'FS=', 'ReadPosRankSum=', 'SOR=']
    # split out all necessary fields from info columns
    parts = dict(part.split('=') for part in info.split(';') if any(field in part for field in fields))
    return parts


def Get_Calls_Info(vcf):
    """Splits and generates necessary info from the vcf calls field

    Parameters
    ----------
    vcf : DataFrame
        Vcf of interest

    Returns
    -------
    DataFrame
        A df of all the fields extracted from the vcf
    """

    call_fields = vcf['CALLS'].str.split(":", expand=True)
    call_fields.columns = vcf['FORMAT'][0].split(":")
    GTS = pd.get_dummies(call_fields['GT']) # Binary rep of hom vs het calls
    AD = call_fields['AD'].str.split(',', expand=True)
    AD.columns = ['RefD', 'AltD', 'AltAltD']
    AD = AD.drop('AltAltD', axis=1)
    AD['RefD'] = pd.to_numeric(AD['RefD']) # Reference allele depth
    AD['AltD'] = pd.to_numeric(AD['AltD']) # Alternate allele depth
    AD['RDper'] = (AD['RefD']/(AD['RefD'] + AD['AltD'])) # percentage ref
    AD['ADrat'] = (AD['AltD']/(AD['RefD'] + .1)) # Ratio of alt to ref
    calls = pd.concat([GTS['0/1'], AD], axis=1)
    return calls


def Check_Type(poly):
    """Ensures input for polymorphism type is an accepted value

    Parameters
    ----------
    poly : str
        A string representing the type of polymorphism. Either 'SNP' or 'INDEL'
    """

    if poly != 'SNP' and poly != 'INDEL':
        raise ValueError('--type takes only values SNP or INDEL')

    return


def Check_VCF(vcf_path):
    """Confirms that input VCF has only one calls column

    Parameters
    ----------
    vcf_path : str
        location of vcf file of interest
    """
    def assert_vcf(vcf):
        reg = re.compile('^#CHROM')
        ms_msg = '''Input VCF has multiple samples.
                    ExtremeVariantFilter only supports single sample VCFs.'''
        no_format_msg = '''Input VCF is missing required FORMAT field'''
        for line in vcf:
            if reg.match(line):
                assert len(line.split('\t')) == 10, ms_msg
                assert "FORMAT" in line.split('\t'), no_format_msg
                break


    if Is_Gzipped(vcf_path):
        with gzip.open(vcf_path, 'r') as vcf:
            assert_vcf(vcf)
    else:
        with open(vcf_path, 'r') as vcf:
            assert_vcf(vcf)


def Is_Gzipped(path):
    """Checks whether a file has a gzip extension or not

    Parameters
    ----------
    vcf_path : str
        location of vcf file of interest
    """
    if path.split('/').pop().split('.')[-1:][0] == 'gz':
        return True
    else:
        return False


# Apply Filter Functions


def Get_Header(vcf_path):
    """Splits header from input vcf and returns it

    Parameters
    ----------
    vcf_path : str
        location of vcf file of interest

    Returns
    -------
    list
        a list representation of the header with added FILTER info
    """

    def read_header(vcf):
        header = []
        newline = vcf.readline()
        while newline.startswith('#'):
            if 'FILTER' in newline or "CHROM\tPOS" in newline and filter_written == False:
                header.append('##FILTER=<ID=EVF_SNP,Description="Likely FP SNP as determined by loaded model">\n')
                header.append('##FILTER=<ID=EVF_IND,Description="Likely FP InDel as determined by loaded model">\n')
                filter_written = True
            header.append(newline)
            newline = vcf.readline()


    if Is_Gzipped(vcf_path):
        with gzip.open(vcf_path, 'rb') as vcf:
            header = read_header(vcf)
    else:
        with open(vcf_path, 'r') as vcf:
            header = read_header(vcf)
    return header


def Check_SNP(vcf):
    """Checks whether input variant is a SNP and returns either 0 or 1

    Parameters
    ----------
    vcf : DataFrame
        vcf file with variants of interest

    Returns
    -------
    int
        an integer representation of whether or not a variant is a SNP
    """

    if "," in vcf['ALT']: # check for multiple alt fields
        if len(vcf['REF']) == 1 and (len(vcf['ALT'].split(',')[0]) == 1 or \
                                     len(vcf['ALT'].split(',')[1]) == 1):
                                     # Determine whether any alt
                                     # fields are indels
            return 1
        else:
            return 0
    elif len(vcf['REF']) == 1 and len(vcf['ALT']) == 1:
        return 1
    else:
        return 0


def Predict_Var(vcf, snp_mdl, ind_mdl):
    """Predicts the truth value of a variant

    Parameters
    ----------
    vcf : DataFrame
        vcf file with variants of interest
    snp_mdl : XGBClassifier
        classifier object, trained on SNPs
    ind_mdl : XGBClassifier
        classifier object, trained on InDels

    Returns
    -------
    int
        an integer representation of whether or not a variant is predicted
        to be real
    """

    params = vcf.iloc[0:11].values
    if vcf['Is_SNP'] == 1:
        return int(snp_mdl.predict(params[None, :]))
    else:
        return int(ind_mdl.predict(params[None, :]))


def Add_Filter(vcf):
    """Adds filter info based on whether a variant is a snp and its predicted
       truth value

    Parameters
    ----------
    vcf : DataFrame
        vcf file with variants of interest

    Returns
    -------
    str
        'EVF_SNP' if the variant is a SNP and predicted False
        'EVF_IND' if the variant is an InDel and predicted False
        '.' for all others
    """

    if vcf['Is_SNP'] == 1 and vcf['Predict'] == 0:
        return "EVF_SNP"
    elif vcf['Is_SNP'] == 0 and vcf['Predict'] == 0:
        return "EVF_IND"
    else:
        return "."


def Get_Name(path):
    """Get the root file name from input path and append '.filter.vcf'

    Parameters
    ----------
    path : str
        location of vcf file of interest

    Returns
    -------
    str
        root file name of path with .filter.vcf appended
    """

    filename = path.split('/').pop()
    if Is_Gzipped(path):
        basename = '.'.join(filename.split('.')[0:-2])
    else:
        basename = '.'.join(filename.split('.')[0:-1])
    outname = basename + '.filter.vcf'
    return outname


def Is_Gzipped(path):
    ""
    if path.split('/').pop().split('.')[-1:][0] == 'gz':
        return True
    else:
        return False


def Write_VCF(table, header, outname):
    """Write out vcf with filter info added

    Parameters
    ----------
    table : DataFrame
        final dataframe with updated filter info
    header : list
        stored header from original vcf with added FILTER lines
    outname : str
        output name for new vcf
    """

    with open(outname, 'w') as vcf_out:
        for line in header:
            vcf_out.write("%s" % line)
    table.to_csv(outname, sep = '\t', header = False, mode = 'a', index = False)


# Train Model Functions


def Make_Table(vcf, label):
    """Creates a table for model training with a label for its truth value

    Parameters
    ----------
    vcf : str
        file location of vcf of interest
    label : int
        1 for true positive vcf or 0 for false positive vcfs

    Returns
    -------
    DataFrame
        a DataFrame representation of the original vcf with an added truth label
    """

    var_vcf = Open_VCF(vcf)
    info_fields = pd.DataFrame(list(var_vcf['INFO'].apply(Split_Info))).fillna(0.)
    info_fields = info_fields[['QD', 'MQ', 'FS', 'MQRankSum', 'ReadPosRankSum', 'SOR']]
    calls = Get_Calls_Info(var_vcf)
    info_fields = pd.concat([info_fields, calls], axis=1)
    info_fields['label'] = label
    return info_fields


def Get_Training_Table(tp_vcf, fp_vcf):
    """Creates a table for model training with both true and false variants

    Parameters
    ----------
    tp_vcf : str
        file location of a true positive vcf from VCFeval
    fp_vcf : str
        file location of a false positive vcf from VCFeval

    Returns
    -------
    DataFrame
        a DataFrame representation of the input vcfs with truth labels
    """

    tp_snp_vcf = Make_Table(tp_vcf, 1)
    fp_snp_vcf = Make_Table(fp_vcf, 0)
    full_vcf = tp_snp_vcf.append(fp_snp_vcf)
    array = full_vcf.values
    X = array[:,0:11]
    Y = array[:,11]
    return X, Y.astype(int)


def Get_Training_Tables(tp_fp_tup):
    """Creates a table for model training with both true and false variants

    Parameters
    ----------
    tp_fp_tup : str
        a tuple with a true and false positive vcf


    Returns
    -------
    DataFrame
        a DataFrame representation of the input vcfs with truth labels
    """

    tp_vcf, fp_vcf = tp_fp_tup
    Check_VCF(tp_vcf)
    Check_VCF(fp_vcf)
    tp_snp_vcf = Make_Table(tp_vcf, 1)
    fp_snp_vcf = Make_Table(fp_vcf, 0)
    full_vcf = tp_snp_vcf.append(fp_snp_vcf)
    array = full_vcf.values
    X = array[:,0:11]
    Y = array[:,11]
    return X, Y.astype(int)


def Build_Model(poly, njobs):
    """Builds an XGBClassifier model with parameters specific to SNP or InDel
       Training vcfs

    Parameters
    ----------
    poly : str
        type of polymorphsim for training. 'SNP' or 'INDEL'
    njobs : int
        number of threads for model to use during training

    Returns
    -------
    XGBClassifier
        an XGBClassifier object with parameters tuned for 'SNP' or 'INDEL'
    """

    seed = 7
    if poly == "SNP":
        model = ("XGBoost('gbtree', 0.3, 6, 600)",
                 XGBClassifier(n_estimators=600,
                 learning_rate=0.3, max_depth=6, random_state=seed,
                 algorithm='gbtree', objective="binary:logistic",
                 nthread=njobs))
    elif poly == "INDEL":
        model = ("XGBoost('gbtree', 0.3, 6, 1000)",
                 XGBClassifier(n_estimators=1000,
                 learning_rate=0.3, max_depth=6, random_state=seed,
                 algorithm='gbtree', objective="binary:logistic",
                 nthread=njobs))
    else:
        raise ValueError('--type takes only values SNP or INDEL')

    return model


def Check_VCF_Paths(tp_vcf_path, fp_vcf_path):
    """Checks to make sure there are equal numbers of tp and fp vcfs

    Parameters
    ----------
    tp_vcf_path : str
        comma seperated list of TP vcf paths
    fp_vcf_path : str
        comma seperated list of FP vcf paths

    Returns
    -------
    list
        a list of tuples containing paired vcf paths
    """

    if "," in tp_vcf_path or "," in fp_vcf_path:
        if len(tp_vcf_path.split(",")) != len(fp_vcf_path.split(",")):
            raise ValueError('Unequal number of True and False VCFs supplied')

    return zip(tp_vcf_path.split(","), fp_vcf_path.split(","))
