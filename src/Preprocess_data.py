# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 17:05:22 2019

@author: rmagni
"""
import pandas as pd, numpy as np, time, os, logging

def redundancy(mzml):
    a = pd.DataFrame(mzml.groupby(["Modified Sequence", "Charge"])["Retention time"].agg(["min", "max", "count"])).reset_index()
    a.columns = ['Modified Sequence', 'Charge', 'min', 'max', 'Redundancy']
    a["Resolution range"] = a["max"] - a["min"]
    a = a[["Modified Sequence", "Charge", "Redundancy", "Resolution range"]]
    return a


def pre_proccess_Data(filein1, filein2, label, isoname, start):
    fileout = filein1.replace(".mzML.parse.tsv", ".QCData.tsv")
    mzml = pd.read_csv(filein1, sep="\t")
    mzml1 = mzml.loc[mzml["ms level"] == 1]
    if label != 4:
        mzml["Reporter count"] = np.where(mzml["ms level"] == 1, np.nan, mzml[isoname].count(axis=1))
        mzml["Median reporter int"] = mzml[isoname].median(axis=1)
        mzml["Median reporter int[log10]"] = np.log10(mzml["Median reporter int"])
        types = {'Scan Number': 'uint32', 'Retention time': 'float32',
         'ms level': 'uint8', 'Scan header': 'float32', 'Injection time': 'float32',
         'Precursor Scan': 'float32', 'Precursor m/z': 'float32',
         'Charge': 'float32', 'BPC': 'float32', 'TIC': 'float32',
         'm/z count': 'uint32', 'Median Intensity': 'float32',
         'Reporter count': 'float32', 'Median reporter int': 'float32',
         'Median reporter int[log10]': 'float32'}
        columns = ['Scan Number', 'Retention time', 'ms level', 'Scan header', 'Injection time',
         'Precursor Scan', 'Precursor m/z',
         'TopN',
         'Intensity jumps', 'Intensity', 'Intensity[log10]', 'Charge',
         'BPC', 'TIC', 'm/z count', 'Median Intensity',
         'Median Intensity[log10]',
         'Reporter count', 'Median reporter int', 'Median reporter int[log10]',
         'Scan', 'Peptide Sequence',
         'Modified Sequence',
         'Protein', 'Matched fragment ions', 'Total possible number of matched theoretical fragment ions',
         'Missed cleavages',
         'Hyperscore', 'Delta Mass [ppm]']
    else:
        types = {
         'Scan Number': 'uint32', 'Retention time': 'float32',
         'ms level': 'uint8', 'Scan header': 'float32', 'Injection time': 'float32',
         'Precursor Scan': 'float32', 'Precursor m/z': 'float32',
         'Charge': 'float32', 'BPC': 'float32', 'TIC': 'float32',
         'm/z count': 'uint32', 'Median Intensity': 'float32'}
        columns = ['Scan Number', 'Retention time', 'ms level', 'Scan header', 'Injection time',
         'Precursor Scan', 'Precursor m/z', 'TopN',
         'Intensity jumps',
         'Intensity', 'Intensity[log10]', 'Charge', 'BPC', 'TIC',
         'm/z count', 'Median Intensity',
         'Median Intensity[log10]',
         'Scan', 'Peptide Sequence', 'Modified Sequence', 'Protein',
         'Matched fragment ions',
         'Total possible number of matched theoretical fragment ions',
         'Missed cleavages', 'Hyperscore', 'Delta Mass [ppm]']
    mzml["Median Intensity[log10]"] = np.log10(mzml["Median Intensity"])
    mzml["Scan header"] = mzml1.groupby("Scan header").ngroup()
    mzml = mzml.astype(types)
    ls = mzml["Scan Number"].max()
    lsc = mzml["Scan header"].max()
    mzml["Intensity"] = np.where(mzml["ms level"] == 1, mzml["BPC"], mzml["TIC"])
    mzml["Intensity[log10]"] = np.log10(mzml["Intensity"])
    mzml2 = mzml1.loc[mzml["Scan header"] == lsc]
    mzml["TopN"] = (-mzml2["Scan Number"].diff(periods=(-1)) - 1 - lsc).fillna(ls - mzml2["Scan Number"].max())
    mzml["Intensity jumps"] = mzml1["TIC"] / mzml1["TIC"].shift(1)
    mzml = mzml.merge((pd.read_csv(filein2, "\t")), left_on="Scan Number", right_on="Scan", how="left")
    mzml = mzml[columns]
    mzml = mzml.astype({'Scan': '"float32"', 'Matched fragment ions': '"float32"', 'Total possible number of matched theoretical fragment ions': '"float32"',
     'Missed cleavages': '"float32"', 'Hyperscore': '"float32"', 'Delta Mass [ppm]': '"float32"'})
    mzml["ID efficiency [%]"] = (mzml1["Scan Number"] - 1).astype("float32")
    mzml["ID efficiency [%]"] = mzml["ID efficiency [%]"].replace(to_replace=(np.nan), method="ffill").astype("float32")
    mzml["ID efficiency [%]"] = (mzml.groupby("ID efficiency [%]")["Modified Sequence"].count() / (mzml.groupby("ID efficiency [%]")["Scan Number"].count() - 1) * 100).astype("float32")
    mzml["Matched ions [%]"] = (mzml["Matched fragment ions"] / mzml["Total possible number of matched theoretical fragment ions"] * 100).astype("float32")
    mzml["Acummulated redundancy"] = (mzml.groupby(["Modified Sequence", "Charge"])["Scan"].cumcount() + 1).astype("float32")
    mzml = mzml.merge((redundancy(mzml)), how="left", left_on=["Modified Sequence", "Charge"], right_on=["Modified Sequence", "Charge"])
    mzml["Redundancy"] = np.where(mzml["Acummulated redundancy"] == 1, mzml["Redundancy"], np.nan).astype("float32")
    mzml["Resolution range"] = np.where(mzml["Acummulated redundancy"] == 1, mzml["Resolution range"], np.nan).astype("float32")
    mzml["Resolution range"] = np.where(mzml["Redundancy"] == 1, np.nan, mzml["Resolution range"]).astype("float32")
    mzml["Redundancy"] = np.where(mzml["Redundancy"] == 1, np.nan, mzml["Redundancy"]).astype("float32")
    mzml.to_csv(fileout, sep="\t", index=False)
    end = time.time()
    timer = divmod(end - start, 60)
    etxt = os.path.basename(fileout)
    etxt1 = "Generated in " + str(timer[0]) + " Minutes and " + str(round(timer[1], 4)) + " Seconds"
    logging.info("\n" + "||".center(80, "|") + etxt.center(80, " ") + etxt1.center(80, " ") + "||".center(80, "|"))
    return mzml