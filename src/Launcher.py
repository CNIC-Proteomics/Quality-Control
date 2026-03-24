# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 17:06:24 2019

@author: rmagni
"""
import multiprocessing, logging, argparse, os, sys, glob, concurrent.futures, subprocess, time
import numpy as np
from pathlib import Path
from itertools import repeat
# Custom modules
import iso_quan_and_correction as iqc
import mzml_parser_completo as mpc
import pratiomsfragger as prms
import qc_report as qc
# Ignore numpy warnings
# np.warnings.filterwarnings("ignore")


def run_mzml_parse(filesinraw, folderoutraw, ThermoRawFileParserpath):
    subprocess.run([ThermoRawFileParserpath, filesinraw, folderoutraw, '-f=2', '-m=0'],
                   check=True, shell=True)


def launcher_task(filesin, folders, inpathraw, log):
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

def main(args):
    logging.info("STARTING QUALITY CONTROL WORKFLOW")
    ROOT_DIR = Path(__file__).parent.parent
    multiprocessing.freeze_support()
    with open(args.params, "r") as f:
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
    msfraggerpath = os.path.join(ROOT_DIR, "lib", "MSFragger", "MSFragger.jar")
    if not os.path.isfile(msfraggerpath):
        logging.error("MSFragger not found in " + str(msfraggerpath))
        sys.exit()
    ThermoRawFileParserpath = os.path.join(ROOT_DIR, "lib", "ThermoRawFileParser", "ThermoRawFileParser.exe")
    if not os.path.isfile(ThermoRawFileParserpath):
        logging.error("ThermoRawFileParser not found in " + str(ThermoRawFileParserpath))
        sys.exit()
    log = np.array([[raw_parser, msfragger, pratio, mzml_parser, QC_report], [raw_parser, msfragger, pratio, mzml_parser, QC_report]])
    fileslist, log, folders = launcher_task(args.filesin, args.folders, args.raw, log)
    start = time.time()
    if log[1][0] == 1:
        filesinraw = ["-i=" + file for file in fileslist[0]]
        folderoutraw = "-o=" + folders
        logging.info("Starting ThermoRawFileParser")
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            futures = executor.map(run_mzml_parse,
                                   filesinraw,
                                   repeat(folderoutraw),
                                   repeat(ThermoRawFileParserpath))
            try:
                for _ in futures:
                    pass
            except:
                logging.error("Error running ThermoRawFileParser")
                executor.shutdown(cancel_futures=True)
                sys.exit()
        end = time.time()
        timer = divmod(end - start, 60)
        logging.info("Finished ThermoRawFileParser in " + str(timer[0]) + " Minutes and " + str(round(timer[1], 4)) + " Seconds")
    if log[1][1] == 1:
        try:
            subprocess.run((["java", "-jar", msfraggerpath, args.params] + fileslist[1]), check=True, shell=True)
        except:
            logging.error("Error running MSFragger")
            sys.exit()
    if log[1][2] == 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            futures = executor.map(prms.Pratio, fileslist[2], repeat(decoy_prefix), repeat(deltaMassThreshold), repeat(FDRlvl), repeat(JumpsAreas))
            for _ in futures:
                logging.info(_)
    if log[1][3] == 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            futures =  executor.map(mpc.parser, fileslist[3], repeat(isobaric_labeling), repeat(isotag), repeat(isoname), repeat(isocorrm), repeat(isobaric_labeling_isotopic_correction), repeat(intensity_noise_level1), repeat(intensity_noise_level2), repeat(mz_intensity_array))
            for _ in futures:
                logging.info(_)
    if log[1][4] == 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            futures = executor.map(qc.create_report, fileslist[5], fileslist[4], repeat(isobaric_labeling), repeat(isoname))
            for _ in futures:
                logging.info(_)
    end = time.time()
    timer = divmod(end - start, 60)
    logging.info("FINISHED QUALITY CONTROL WORKFLOW IN " + str(timer[0]) + " MINUTES AND " + str(round(timer[1], 4)) + " SECONDS")
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Quality-Control',
        epilog='''
        Example:
            python Quality-Control.py -r raw -p params

        ''')
    multiprocessing.freeze_support()
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--raw",    required=True, help="raw file or raw folder", type=Path)
    parser.add_argument("-p", "--params", required=True, help="parameters file",        type=Path)
    args = parser.parse_args()
    
    args.raw = str(args.raw)
    args.params = str(args.params)
    if args.raw.endswith(".raw") == True:
        filesin = [args.raw]
        folders = os.path.join(os.path.dirname(args.raw), "Data")
    else:
        filesin = glob.glob(os.path.join(args.raw, "*.raw"))
        folders = os.path.join(args.raw, "Data")
    if not os.path.exists(folders):
        os.makedirs(folders)
    args.filesin = filesin
    args.folders = folders
    
    log_file = os.path.join(args.folders + '/Quality-Control.log')
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        handlers=[logging.FileHandler(log_file),
                                  logging.StreamHandler()])
    
    if len(args.filesin) < 1:
        logging.warning("No RAW files found in " + str(args.raw))
        sys.exit()
    
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')