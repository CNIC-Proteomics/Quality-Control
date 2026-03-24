# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 16:12:20 2026

@author: alaguillog
"""

import subprocess

def run_mzml_parse(filesinraw, folderoutraw, ThermoRawFileParserpath):
    subprocess.run([ThermoRawFileParserpath, filesinraw, folderoutraw, '-f=2', '-m=0'],
                   check=True, shell=True)