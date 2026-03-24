# -*- coding: utf-8 -*-
"""
Created on Mon May 27 15:33:42 2019

@author: rmagni
"""
import time, itertools, pandas as pd, numpy as np, os

h_mass = 1.00727647

def Jumps(df, JumpsAreas):
    s1 = df["Mass difference"] / (df["Neutral mass of peptide (including any variable modifications) (Da)"] + h_mass) * 1000000.0
    s2 = (df["Mass difference"] - 1.0033) / (df["Neutral mass of peptide (including any variable modifications) (Da)"] + h_mass) * 1000000.0
    s3 = (df["Mass difference"] + 1.0033) / (df["Neutral mass of peptide (including any variable modifications) (Da)"] + h_mass) * 1000000.0
    s4 = (df["Mass difference"] - 2.0066) / (df["Neutral mass of peptide (including any variable modifications) (Da)"] + h_mass) * 1000000.0
    s5 = (df["Mass difference"] + 2.0066) / (df["Neutral mass of peptide (including any variable modifications) (Da)"] + h_mass) * 1000000.0
    df1 = pd.DataFrame({'1': s1, '2': s2, '3': s3, '4': s4, '5': s5})
    if JumpsAreas == "1":
        Jumps = ["1"]
    elif JumpsAreas == "3":
        Jumps = ["1", "2", "3"]
    else:
        if JumpsAreas == "5":
            Jumps = ['1', '2', '3', '4', '5']
    Deltamass = df1[Jumps].abs().min(axis=1)
    return Deltamass


def targetdecoy(df, decoy_prefix):
    z = list(df["Protein"].str.split(";"))
    p = [all(decoy_prefix in item for item in i) for i in z]
    p = list(map(int, p))
    return p


def SequenceMod(df):
    s = list(df["Variable modifications detected"].fillna("0A").replace({'N_term(.+?),':"",  '([A-Z])':";"}, regex=True).str.split(","))
    s = [[a.split(";") for a in i] for i in s]
    sn = [[a for a, b in i] for i in s]
    si = [[b for a, b in i] for i in s]
    sn = [list(filter(None, i)) for i in sn]
    sn = [[int(i) for i in j] for j in sn]
    [[a.insert(0, 0)] for a in sn]
    l = list(df["Peptide Sequence"])
    f = [[a[i:j] for i, j in zip(b, b[1:] + [None])] for a, b in zip(l, sn)]
    x = ["".join(list(itertools.chain.from_iterable(list(itertools.zip_longest(i, j, fillvalue=""))))) for i, j in list(zip(f, si))]
    return x


def FdrXc(df, FDRlvl):
    df = df.sort_values(by="Hyperscore", ascending=False)
    df["rank"] = df.groupby("T_D").cumcount() + 1
    df["rank_T"] = np.where(df["T_D"] == 0, df["rank"], 0)
    df["rank_T"] = df["rank_T"].replace(0, np.nan).ffill()
    df["rank_D"] = np.where(df["T_D"] == 1, df["rank"], 0)
    df["rank_D"] = df["rank_D"].replace(0, np.nan).ffill()
    df["rank"] = df.index + 1
    df["FdrXc"] = df["rank_D"] / df["rank_T"]
    df = df[df["FdrXc"] <= FDRlvl]
    df = df[df["T_D"] == 0]
    return df


def Pratio(filein, decoy_prefix, deltaMassThreshold, FDRlvl, JumpsAreas):
    fileout = filein.replace(".result.tsv", ".filtered.result.tsv")
    # stxt = "Processing " + os.path.basename(filein)
    # logging.info("\n" + "*".center(80, "*") + stxt.center(80, " ") + "*".center(80, "*"))
    start = time.time()
    df = pd.read_csv(filein, sep="\t")
    df.columns = ['Scan',
     'Precursor neutral mass (Da)',
     'Retention time (minutes)',
     'Precursor charge',
     'Hit rank',
     'Peptide Sequence',
     'Upstream Amino Acid',
     'Downstream Amino Acid',
     'Protein',
     'Matched fragment ions',
     'Total possible number of matched theoretical fragment ions',
     'Neutral mass of peptide (including any variable modifications) (Da)',
     'Mass difference',
     'Number of tryptic termini',
     'Missed cleavages',
     'Variable modifications detected',
     'Hyperscore',
     'Next score',
     'Intercept of expectation model (expectation in log space)',
     'Slope of expectation model (expectation in log space)',
     'best_locs',
     'score_without_delta_mass',
     'best_score_with_delta_mass',
     'second_best_score_with_delta_mass']
    df["Delta Mass [ppm]"] = Jumps(df, JumpsAreas)
    df = df[df["Delta Mass [ppm]"].abs() <= deltaMassThreshold]
    df["T_D"] = targetdecoy(df, decoy_prefix)
    df = FdrXc(df, FDRlvl)
    df["Modified Sequence"] = SequenceMod(df)
    df.to_csv(fileout, sep="\t")
    end = time.time()
    timer = divmod(end - start, 60)
    return str(os.path.basename(fileout)) + " generated in " + str(timer[0]) + " Minutes and " + str(round(timer[1], 4)) + " Seconds"