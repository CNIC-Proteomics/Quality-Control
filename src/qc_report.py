# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 14:10:45 2019

@author: rmagni
"""
from matplotlib.backends.backend_pdf import PdfPages
import grc, Preprocess_data as prd, time, os
import logging

def create_pdf(mzml, fileout, label):
    mzml = mzml
    pp = PdfPages(fileout)
    grc.main_report(mzml, pp, label)
    grc.sampling_report(mzml, pp)
    grc.hist_report(mzml, 1, 0, "Intensity", pp)
    grc.hist_report(mzml, 1, 0, "Intensity[log10]", pp)
    grc.hist_report(mzml, 1, 0, "Median Intensity", pp)
    grc.hist_report(mzml, 1, 0, "Median Intensity[log10]", pp)
    grc.hist_report(mzml, 1, 0, "Injection time", pp)
    grc.hist_report(mzml, 1, 0, "m/z count", pp)
    grc.hist_report(mzml, 2, 1, "Intensity", pp)
    grc.hist_report(mzml, 2, 1, "Intensity[log10]", pp)
    grc.hist_report(mzml, 2, 1, "Median Intensity", pp)
    grc.hist_report(mzml, 2, 1, "Median Intensity[log10]", pp)
    grc.hist_report(mzml, 2, 1, "Injection time", pp)
    grc.hist_report(mzml, 2, 1, "m/z count", pp)
    grc.bar_report(mzml, 1, 0, "TopN", pp)
    grc.bar_report(mzml, 2, 1, "Charge", pp)
    grc.hist_report(mzml, 2, 1, "Precursor m/z", pp)
    if label != 4:
        grc.bar_report(mzml, 2, 1, "Reporter count", pp)
        grc.hist_report(mzml, 2, 1, "Median reporter int", pp)
        grc.hist_report(mzml, 2, 1, "Median reporter int[log10]", pp)
    grc.hist_report(mzml, 1, 0, "ID efficiency [%]", pp)
    grc.hist_report(mzml, 2, 0, "Delta Mass [ppm]", pp)
    grc.hist_report(mzml, 2, 0, "Hyperscore", pp)
    grc.hist_report(mzml, 2, 0, "Matched ions [%]", pp)
    grc.bar_report(mzml, 2, 0, "Missed cleavages", pp)
    pp.close()


def create_report(filein1, filein2, label, isoname):
    stxt = "Processing " + os.path.basename(filein1)
    stxt1 = os.path.basename(filein2)
    logging.info("\n" + "**".center(80, "*") + stxt.center(80, " ") + "and".center(80, " ") + stxt1.center(80, " ") + "**".center(80, "*"))
    start = time.time()
    fileout = filein1.replace(".mzML.parse.tsv", ".QCreport.pdf")
    mzml = prd.pre_proccess_Data(filein1, filein2, label, isoname, start)
    create_pdf(mzml, fileout, label)
    end = time.time()
    timer = divmod(end - start, 60)
    etxt = os.path.basename(fileout)
    etxt1 = "Generated in " + str(timer[0]) + " Minutes and " + str(round(timer[1], 4)) + " Seconds"
    logging.info("\n" + "||".center(80, "|") + etxt.center(80, " ") + etxt1.center(80, " ") + "||".center(80, "|"))