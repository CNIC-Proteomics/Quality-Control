# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 15:08:13 2019

@author: rmagni
"""
import numpy as np, pandas as pd

def isobaric_labelling(label):
    if label == 0:
        p = [[126.127726, 127.124761, 127.131081, 128.128116, 128.134436, 129.131471,
          129.13779, 130.134825, 130.141145, 131.13818],
         [
          '126', '127N', '127C', '128N', '128C', '129N',
          '129C', '130N', '130C', '131']]
    elif label == 1:
        p = [[126.127726, 127.124761, 128.134436, 129.131471, 130.141145, 131.13818], ['126', '127', '128', '129', '130', '131']]
    elif label == 2:
        p = [[113.1078, 114.1112, 115.1082, 116.1116, 117.1149, 118.112, 119.1153,
          121.122], ['113', '114', '115', '116', '117', '118', '119', '121']]
    elif label == 3:
        p = [[114.1112, 115.1083, 116.1116, 117.115], ["114", "115", "116", "117"]]
    else:
        if label == 4:
            p = [[None], [None]]
    return p


def correcmatrix(m, cor):
    if cor == 1:
        m = pd.read_excel(m)
        m.columns = ['Mass Tag', 'Reporter ion', '-2', '-1', 'Monoisotopic', '+1', '+2']
        m = m[["-2", "-1", "+1", "+2"]].values.T
        isocorrm = np.array([[100, 0, m[1][2], 0, m[0][4], 0, 0, 0, 0, 0],
         [
          0, 100, 0, m[1][3], 0, m[0][5], 0, 0, 0, 0],
         [
          m[2][0], 0, 100, 0, m[1][4], 0, m[0][6], 0, 0, 0],
         [
          0, m[2][1], 0, 100, 0, m[1][5], 0, m[0][7], 0, 0],
         [
          m[3][0], 0, m[2][2], 0, 100, 0, m[1][6], 0, m[0][8], 0],
         [
          0, m[3][1], 0, m[2][3], 0, 100, 0, m[1][7], 0, m[0][9]],
         [
          0, 0, m[3][2], 0, m[2][4], 0, 100, 0, m[1][8], 0],
         [
          0, 0, 0, m[3][3], 0, m[2][5], 0, 100, 0, m[1][9]],
         [
          0, 0, 0, 0, m[3][4], 0, m[2][6], 0, 100, 0],
         [
          0, 0, 0, 0, 0, m[3][5], 0, m[2][7], 0, 100]]) / 100
    else:
        isocorrm = None
    return isocorrm


def monoisocorrec(b1, isocorrm):
    b1 = np.array(b1)
    bid = np.array([i != "" for i in b1])
    c1 = [i for i, x in enumerate(bid) if not x]
    b1 = np.linalg.solve(isocorrm[np.ix_(bid, bid)], b1[bid].astype("float32"))
    for i in c1:
        b1 = np.insert(b1, i, np.nan)

    return b1


def get_quant(a, b, label, isotag, isocorrm):
    ppm = 5
    p1 = []
    b1 = []
    po = np.zeros((len(a)), dtype=bool)
    for i, j in enumerate(isotag):
        if np.any(np.logical_and(np.greater_equal(a, j - j * ppm * 1e-06), np.less_equal(a, j + j * ppm * 1e-06))) == True:
            po += np.logical_and(np.greater_equal(a, j - j * ppm * 1e-06), np.less_equal(a, j + j * ppm * 1e-06))
            p1.append(i)

    a = a[po]
    b = b[po]
    p1 = np.array(p1)
    for i in np.arange(len(isotag)):
        b2 = b[np.where(p1 == i)]
        if len(b2) == 0:
            b1.append("")
        else:
            b1.append(b2[0])

    if label == 0:
        if isinstance(isocorrm, np.ndarray):
            b1 = monoisocorrec(b1, isocorrm)
    return b1