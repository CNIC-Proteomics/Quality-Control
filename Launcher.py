# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 17:06:24 2019

@author: rmagni
"""
import multiprocessing
from pathlib import Path
from itertools import repeat
import argparse, os, sys, glob, numpy as np, concurrent.futures, iso_quan_and_correction as iqc, mzml_parser_completo as mpc, pratiomsfragger as prms, qc_report as qc, subprocess, time
np.warnings.filterwarnings("ignore")

def run_mzml_parse(filesinraw, folderoutraw, ThermoRawFileParserpath):
    subprocess.run([ThermoRawFileParserpath, filesinraw, folderoutraw, '-f=2', '-m=0'], shell=True)


def launcher_task(inpathraw, log):
    if inpathraw.endswith(".raw") == True:
        filesin = [inpathraw]
        folders = os.path.join(os.path.dirname(inpathraw), "Data")
    else:
        filesin = glob.glob(os.path.join(inpathraw, "*.raw"))
        folders = os.path.join(inpathraw, "Data")
    if not os.path.exists(folders):
        os.makedirs(folders)
    filesmzML = [os.path.join(os.path.dirname(i), folders, os.path.basename(i).replace(".raw", ".mzML")) for i in filesin]
    filesresult = [i.replace(".mzML", ".result.tsv") for i in filesmzML]
    filesfilteredresult = [i.replace(".mzML", ".filtered.result.tsv") for i in filesmzML]
    filesmzMLparse = [i + ".parse.tsv" for i in filesmzML]
    fileslist = [filesin, filesmzML, filesresult, filesmzML, filesfilteredresult, filesmzMLparse]
    dictfile = {'a':[],  'b':[],  'c':[],  'd':[],  'e':[]}
    if log[(1, 4)] == 1:
        dictfile["a"] = [str(i) for i in np.array(filesresult)[np.array([not os.path.isfile(i) for i in filesfilteredresult])]]
        dictfile["b"] = [str(i) for i in np.array(filesmzML)[np.array([not os.path.isfile(i) for i in filesmzMLparse])]]
        if not dictfile["a"] == []:
            log[(1, 2)] = 1
        if not dictfile["b"] == []:
            log[(1, 3)] = 1
    if log[(1, 2)] == 1:
        dictfile["c"] = [str(i) for i in np.array(filesmzML)[np.array([not os.path.isfile(i) for i in filesresult])]]
        if not dictfile["c"] == []:
            log[(1, 1)] = 1
    if log[(1, 1)] == 1:
        dictfile["d"] = [str(i) for i in np.array(filesin)[np.array([not os.path.isfile(i) for i in filesmzML])]]
        if not dictfile["d"] == []:
            log[(1, 0)] = 1
    if log[(1, 3)] == 1:
        dictfile["e"] = [str(i) for i in np.array(filesin)[np.array([not os.path.isfile(i) for i in filesmzML])]]
        if not dictfile["e"] == []:
            log[(1, 0)] = 1
    dictfile["f"] = [str(i) for i in np.unique(dictfile["d"] + dictfile["e"])]
    fileslistcorrection = [
     dictfile["f"], dictfile["c"], dictfile["a"], dictfile["b"]]
    for num, i in enumerate(zip(log[0], log[1])):
        if i[0] == 0 and i[1] == 1:
            fileslist[num] = fileslistcorrection[num]

    return (fileslist, log, folders)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--raw", help="raw file or raw folder", type=Path)
    parser.add_argument("-p", "--params", help="params file", type=Path)
    args = parser.parse_args()
    inpathraw = str(args.raw)
    paramsfile = str(args.params)
    with open(paramsfile, "r") as f:
        x = f.read().splitlines()
    x = [a.replace(" ", "") for a in x]
    x = [a.replace("\t", "") for a in x]
    x = dict(list(filter(None, [list(filter(None, s.split("#")[0].split("="))) for s in x])))
    raw_parser = int(x["raw_parser"])
    msfragger = int(x["msfragger"])
    pratio = int(x["pratio"])
    mzml_parser = int(x["mzml_parser"])
    QC_report = int(x["QC_report"])
    intensity_noise_level1 = int(x["intensity_noise_level(ms1)"])
    intensity_noise_level2 = int(x["intensity_noise_level(ms2)"])
    mz_intensity_array = int(x["m/z_intensity_array"])
    isobaric_labeling = int(x["isobaric_labeling"])
    isobaric_labeling_isotopic_correction = int(x["isobaric_labeling_isotopic_correction"])
    isotopic_Distribution_table = x["isotopic_Distribution_table"]
    decoy_prefix = x["decoy_prefix"]
    deltaMassThreshold = float(x["deltaMassThreshold"])
    FDRlvl = float(x["FDRlvl"])
    JumpsAreas = x["JumpsAreas"]
    num_threads = int(x["num_threads"])
    num_threads = num_threads if num_threads > 0 else None
    isocorrm = iqc.correcmatrix(isotopic_Distribution_table, isobaric_labeling_isotopic_correction)
    isotag, isoname = iqc.isobaric_labelling(isobaric_labeling)
    msfraggerpath = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), "msfragger", "MSFragger.jar")
    ThermoRawFileParserpath = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), "ThermoRawFileParser", "ThermoRawFileParser.exe")
    log = np.array([[raw_parser, msfragger, pratio, mzml_parser, QC_report], [raw_parser, msfragger, pratio, mzml_parser, QC_report]])
    fileslist, log, folders = launcher_task(inpathraw, log)
    stxt = "STARTING QUALITY CONTROL WORKFLOW"
    print("\n" + "##".center(80, "#") + stxt.center(80, " ") + "##".center(80, "#"))
    start = time.time()
    if log[1][0] == 1:
        filesinraw = ["-i=" + file for file in fileslist[0]]
        folderoutraw = "-o=" + folders
        stxt = "Starting ThermoRawFileParser"
        print("\n" + "**".center(80, "*") + stxt.center(80, " ") + "**".center(80, "*"))
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            executor.map(run_mzml_parse, filesinraw, repeat(folderoutraw), repeat(ThermoRawFileParserpath))
        end = time.time()
        timer = divmod(end - start, 60)
        etxt1 = "Finished ThermoRawFileParser in " + str(timer[0]) + " Minutes and " + str(round(timer[1], 4)) + " Seconds"
        print("\n" + "||".center(80, "|") + etxt1.center(80, " ") + "||".center(80, "|"))
    if log[1][1] == 1:
        subprocess.run((["java", "-jar", msfraggerpath, paramsfile] + fileslist[1]), shell=True)
    if log[1][2] == 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            executor.map(prms.Pratio, fileslist[2], repeat(decoy_prefix), repeat(deltaMassThreshold), repeat(FDRlvl), repeat(JumpsAreas))
    if log[1][3] == 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            executor.map(mpc.parser, fileslist[3], repeat(isobaric_labeling), repeat(isotag), repeat(isoname), repeat(isocorrm), repeat(isobaric_labeling_isotopic_correction), repeat(intensity_noise_level1), repeat(intensity_noise_level2), repeat(mz_intensity_array))
    if log[1][4] == 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            executor.map(qc.create_report, fileslist[5], fileslist[4], repeat(isobaric_labeling), repeat(isoname))
    end = time.time()
    timer = divmod(end - start, 60)
    etxt1 = "FINISHED QUALITY CONTROL WORKFLOW IN " + str(timer[0]) + " MINUTES AND " + str(round(timer[1], 4)) + " SECONDS"
    print("\n" + "##".center(80, "#") + etxt1.center(80, " ") + "##".center(80, "#"))