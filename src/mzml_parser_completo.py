# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:10:16 2019

@author: rmagni
"""
import time, pandas as pd
from lxml import etree as ET
import numpy as np, base64, zlib, iso_quan_and_correction as iqc, os
import logging

def fast_iter(context, func, fh, label, isotag, isocorrm, nl1, nl2, array):
    for event, elem in context:
        fh.append(func(elem, label, isotag, isocorrm, nl1, nl2, array))
        elem.clear()
        for ancestor in elem.xpath("ancestor-or-self::*"):
            while ancestor.getprevious() is not None:
                del ancestor.getparent()[0]

    del context


def array_decoder(a, ctype, dtypea):
    a = base64.b64decode(a)
    if ctype == "zlib compression":
        a = np.frombuffer((zlib.decompress(a)), dtype=dtypea)
    else:
        a = np.frombuffer(a, dtype=dtypea)
    return a.astype("float32")


def get_spectrum_values(elt, label, isotag, isocorrm, nl1, nl2, array):
    val = []
    type_array = {'64-bit float':np.float64,
     '32-bit float':np.float32}
    class_array = {'m/z array':"mz",  'intensity array':"i"}
    spec = {'sn':int(elt.attrib["id"].split("scan=")[1].split(" ")[0]),
     'RT':"",
     'mslevel':"",
     'sc':"2",
     'IT':"",
     'psn':"",
     'pmz':"",
     'Charge':"",
     'BPC':"",
     'TIC':"",
     'lmz':"",
     'mi':"",
     'mz':"",
     'i':""}
    for elem in elt.iter():
        if elem.tag == "{http://psi.hupo.org/ms/mzml}cvParam":
            if elem.attrib["name"] == "ms level":
                spec["mslevel"] = int(elem.attrib["value"])
            else:
                if elem.attrib["name"] == "filter string":
                    spec["sc"] = str(elem.attrib["value"])
                else:
                    if elem.attrib["name"] == "base peak intensity":
                        spec["BPC"] = round(float(elem.attrib["value"]), 4)
                    else:
                        if elem.attrib["name"] == "total ion current":
                            spec["TIC"] = round(float(elem.attrib["value"]), 4)
                        if spec["mslevel"] == 2:
                            if elem.attrib["name"] == "selected ion m/z":
                                spec["pmz"] = round(float(elem.attrib["value"]), 4)
                        if elem.attrib["name"] == "scan start time":
                            spec["RT"] = round(float(elem.attrib["value"]), 4)
                    if elem.attrib["name"] == "ion injection time":
                        spec["IT"] = round(float(elem.attrib["value"]), 4)
                if spec["mslevel"] == 2:
                    if elem.attrib["name"] == "charge state":
                        spec["Charge"] = int(elem.attrib["value"])
                if elem.tag == "{http://psi.hupo.org/ms/mzml}precursor":
                    spec["psn"] = int(elem.attrib["spectrumRef"].split("scan=")[1].split(" ")[0])
            if elem.tag == "{http://psi.hupo.org/ms/mzml}binaryDataArray":
                val = [elem[0].attrib["name"], elem[1].attrib["name"], elem[2].attrib["name"]]
                ctype = ["zlib compression" if "zlib compression" in val else ""][0]
                dtypea = [type_array[k] for k in val if k in type_array][0]
                vname = [class_array[k] for k in val if k in class_array][0]
                spec[vname] = elem[3].text
                spec[vname] = array_decoder(spec[vname], ctype, dtypea)

    if spec["mslevel"] == 2:
        if label != 4:
            for i, a in enumerate(iqc.get_quant(spec["mz"], spec["i"], label, isotag, isocorrm)):
                spec[i] = a

    if spec["mslevel"] == 1:
        spec["mz"] = spec["mz"][spec["i"] >= nl1]
        spec["i"] = spec["i"][spec["i"] >= nl1]
    if spec["mslevel"] == 2:
        spec["mz"] = spec["mz"][spec["i"] >= nl2]
        spec["i"] = spec["i"][spec["i"] >= nl2]
    spec["lmz"] = len(spec["mz"])
    spec["mi"] = np.median(spec["i"])
    if array == 0:
        spec["mz"] = ""
        spec["i"] = ""
    return list(spec.values())


def parser(filein, label, isotag, isoname, isocorrm, cor, nl1, nl2, array):
    fileout = filein + ".parse.tsv"
    stxt = "Processing " + os.path.basename(filein)
    logging.info("\n" + "**".center(80, "*") + stxt.center(80, " ") + "**".center(80, "*"))
    start = time.time()
    fh = []
    context = ET.iterparse(filein, tag="{http://psi.hupo.org/ms/mzml}spectrum")
    fast_iter(context, get_spectrum_values, fh, label, isotag, isocorrm, nl1, nl2, array)
    columns = ['Scan Number', 'Retention time', 'ms level', 'Scan header', 'Injection time',
     'Precursor Scan', 'Precursor m/z',
     'Charge', 'BPC', 'TIC',
     'm/z count', 'Median Intensity', 'm/z array', 'Intensity array'] + isoname
    columns = list(filter(None, columns))
    fh = pd.DataFrame(fh, columns=columns)
    if array == 0:
        fh = fh.drop(["m/z array", "Intensity array"], axis=1)
    fh.to_csv(fileout, "\t")
    end = time.time()
    timer = divmod(end - start, 60)
    etxt = os.path.basename(fileout)
    etxt1 = "Generated in " + str(timer[0]) + " Minutes and " + str(round(timer[1], 4)) + " Seconds"
    logging.info("\n" + "||".center(80, "|") + etxt.center(80, " ") + etxt1.center(80, " ") + "||".center(80, "|"))