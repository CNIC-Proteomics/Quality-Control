# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 14:12:44 2019

@author: rmagni
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter, MaxNLocator
from matplotlib.font_manager import FontProperties
plt.style.use("seaborn-v0_8-ticks")

class ScalarFormatterForceFormat(ScalarFormatter):

    def _set_format(self, vmin, vmax):
        if vmax < 1000:
            self.format = "%1.1f"
        else:
            self.format = "%1.2f"


class Labeloffset:

    def __init__(self, ax, label='', axis='y'):
        self.axis = {'y':ax.yaxis,
         'x':ax.xaxis}[axis]
        self.label = label
        ax.callbacks.connect(axis + "lim_changed", self.update)
        ax.figure.canvas.draw()
        self.update(None)

    def update(self, lim):
        fmt = self.axis.get_major_formatter()
        self.axis.offsetText.set_visible(False)
        self.axis.set_label_text((self.label + " " + fmt.get_offset()), fontweight="bold")


def time_atr(mzml):
    rtmin = mzml["Retention time"].min()
    rtmax = mzml["Retention time"].max()
    rt_total = round(mzml["Retention time"].max() - mzml["Retention time"].min(), 4)
    return (rtmin, rtmax, rt_total)


def outQ(data, axis=None):
    a, b = np.percentile(data, q=[10, 90], out=None, overwrite_input=False, interpolation="linear", keepdims=False)
    c = b - a
    out = data[data < b + 1.5 * c]
    return out


def bins_range(nmax, nmin):
    a = ""
    nmax = round(nmax / 2) * 2
    nmin = round(nmin / 2) * 2
    if nmax != int(nmax) or nmin != int(nmin):
        nmax = nmax * 10
        nmin = nmin * 10
        a = 1
    else:
        n = nmax - nmin
        if n < 10:
            n = int(n)
            f = [10, 5, 4, 3, 2, 2, 2, 2, 2]
            factor = factors(n)
            factor = [i * f[n - 1] for i in factors(n)]
            factor = [x for x in factor if x >= 10 if x < 21]
        else:
            factor = factors(n)
            factor = [x for x in factor if x >= 10 if x < 21]
            factor = [x for x in factor if n / x % 2 == 0 or n / x == 1]
            while not factor:
                n = n + 1
                factor = factors(n)
                factor = [x for x in factor if x >= 10 if x < 21]
                factor = [x for x in factor if n / x % 2 == 0 or n / x == 1]

        if a == 1:
            n = n / 10
            nmax = nmax / 10
            nmin = nmin / 10
        if n + nmin > nmax:
            if nmin - (n + nmin - nmax) < 0:
                nmin = 0
            else:
                nmin = nmin - (n + nmin - nmax)
            nmax = n + nmin
        else:
            nmax = n + nmin
    return (
     max(factor), nmax, nmin)


def factors(n):
    factors = []
    for i in range(1, 21):
        if n % i == 0:
            factors.append(i)

    return factors


def quantiles(mzml):
    maxid = mzml["Scan"].idxmax()
    minid = mzml.index[mzml["Scan Number"] == int(mzml["Precursor Scan"][mzml["Scan"] == mzml["Scan"].min()].iloc[0])][0]
    q1 = mzml.index[mzml["Scan Number"] == int(mzml["Precursor Scan"][mzml["Scan"] == int(mzml["Scan"].quantile(0.25, interpolation="lower"))].iloc[0])][0]
    q2 = mzml.index[mzml["Scan Number"] == int(mzml["Precursor Scan"][mzml["Scan"] == int(mzml["Scan"].quantile(0.5, interpolation="lower"))].iloc[0])][0]
    q3 = mzml.index[mzml["Scan Number"] == int(mzml["Precursor Scan"][mzml["Scan"] == int(mzml["Scan"].quantile(0.75, interpolation="lower"))].iloc[0])][0]
    id_range = slice(minid, maxid)
    aq_range = slice(minid, q1 - 1)
    bq_range = slice(q1, q2 - 1)
    cq_range = slice(q2, q3 - 1)
    dq_range = slice(q3, maxid)
    return (
     id_range, aq_range, bq_range, cq_range, dq_range)


def rangescal(mzml):
    id_range, aq_range, bq_range, cq_range, dq_range = quantiles(mzml)
    rangesl = zip([[time_atr(mzml.loc[aq_range])[0], time_atr(mzml.loc[dq_range])[1]],
     [
      time_atr(mzml.loc[aq_range])[0], time_atr(mzml.loc[aq_range])[1]],
     [
      time_atr(mzml.loc[bq_range])[0], time_atr(mzml.loc[bq_range])[1]],
     [
      time_atr(mzml.loc[cq_range])[0], time_atr(mzml.loc[cq_range])[1]],
     [
      time_atr(mzml.loc[dq_range])[0], time_atr(mzml.loc[dq_range])[1]]], [
     0.5, 0, 0, 0, 0], [
     1, 0.5, 0.5, 0.5, 0.5], [
     '#ffffcc', '#c2ccff', '#ffa756', '#b1d27b', '#d9544d'], [
     0.75, 0.25, 0.25, 0.25, 0.25], [
     'ID Range(IDR)', 'R1', 'R2', 'R3', 'R4'])
    quantl = zip([time_atr(mzml.loc[aq_range])[0],
     time_atr(mzml.loc[bq_range])[0],
     time_atr(mzml.loc[cq_range])[0],
     time_atr(mzml.loc[dq_range])[0],
     time_atr(mzml.loc[dq_range])[1]], [
     (0, (2, 4)), (0, (2, 4)), (0, (2, 4)), (0, (2, 4)), (0, (2, 4))], [
     'b', 'orange', 'green', 'r', 'k'], [
     'FirstID', 'ID Q1', 'ID Q2', 'ID Q3', 'LastID'])
    return (
     id_range, aq_range, bq_range, cq_range, dq_range, rangesl, quantl)


def hist_report(mzml, l, p, t, pp):
    id_range, aq_range, bq_range, cq_range, dq_range, rangesl, quantl = rangescal(mzml)
    maxrt = int(mzml["Retention time"].values.max())
    mzmla = mzml.loc[(mzml["ms level"] == l, ["Scan Number", "Retention time", "BPC", "TIC", t, "Scan"])]
    if l == 1:
        mzmlo = mzmla.loc[id_range]
        data = list(zip([mzmla], ["Scans"]))
        datah = list(zip([mzmlo], ["Scans"]))
        cc = ['#ffffcc', '#c2ccff', '#ffa756', '#b1d27b', '#d9544d']
        dmap = ["#0000FF"]
    elif l == 2 and p == 1:
        mzml1 = mzmla[mzmla["Scan"].isna()]
        mzml2 = mzmla.dropna(0, subset=["Scan"])
        data = list(zip([mzml1, mzml2], ["Scans", "PSMs"]))
        datah = list(zip([mzmla, mzml2], ["Scans", "PSMs"]))
        cc = ['#ffffcc', '#cccc99', '#c2ccff', '#8f99cc', '#ffa756', '#cc7423', '#b1d27b',
         '#7e9f48', '#d9544d', '#a6211a']
        dmap = ["#0000FF", "#FF0000"]
    elif l == 2:
        pass
    if p == 0:
        mzmla = mzmla.dropna(0, subset=["Scan"])
        data = list(zip([mzmla], ["PSMs"]))
        datah = list(zip([mzmla], ["PSMs"]))
        cc = ['#ffffcc', '#c2ccff', '#ffa756', '#b1d27b', '#d9544d']
        dmap = ["#FF0000"]
    else:
        c = [
         [
          '#ffffcc', '#c2ccff', '#ffa756', '#b1d27b', '#d9544d'], ['#cccc99', '#8f99cc', '#cc7423', '#7e9f48', '#a6211a']]
        ec = ['#80804d', '#434d80', '#802800', '#325300', '#5a0000']
        TBP = ["", "BPC", "TIC"]
        legendlab = ["Scans Relative Frequency(%)", "PSMs Relative Frequency(%)"]
        ranges = list(zip([id_range, aq_range, bq_range, cq_range, dq_range], ['IDR', 'R1', 'R2', 'R3', 'R4']))
        nmax = round(outQ(datah[0][0][t].dropna()).max(), 1)
        nmin = round(outQ(datah[0][0][t].dropna()).min(), 1)
        bins, nmax, nmin = bins_range(nmax, nmin)
        lgcolor = ["#0000FF", "#FF0000"]
        fig = plt.figure(figsize=(16, 9))
        gs1 = gridspec.GridSpec(12, 1)
        gs1.update(left=0.05, right=0.52, hspace=0.0, top=0.95, bottom=0.4)
        ax1 = fig.add_subplot(gs1[2:8])
        ax11 = fig.add_subplot((gs1[8:12]), sharex=ax1)
        ax12 = fig.add_subplot((gs1[1:2]), sharex=ax1)
        ax13 = fig.add_subplot((gs1[0:1]), sharex=ax1)
        gs2 = gridspec.GridSpec(5, 1)
        gs2.update(left=0.58, right=0.98, top=0.98, bottom=0.08, hspace=0)
        ax21 = fig.add_subplot(gs2[0])
        ax21.xaxis.set_visible(False)
        ax22 = fig.add_subplot((gs2[1]), sharex=ax21)
        ax22.xaxis.set_visible(False)
        ax23 = fig.add_subplot((gs2[2]), sharey=ax22, sharex=ax22)
        ax23.xaxis.set_visible(False)
        ax24 = fig.add_subplot((gs2[3]), sharey=ax22, sharex=ax22)
        ax24.xaxis.set_visible(False)
        ax25 = fig.add_subplot((gs2[4]), sharey=ax22, sharex=ax22)
        ax1.set_xlim(0, maxrt)
        ax1.xaxis.set_visible(False)
        ax12.xaxis.set_visible(False)
        ax12.yaxis.set_visible(False)
        ax13.xaxis.set_visible(False)
        ax13.yaxis.set_visible(False)
        r = [ax21, ax22, ax23, ax24, ax25]
        gs22 = gridspec.GridSpec(1, 1)
        gs22.update(left=0.06, right=0.52, top=0.32, bottom=0.05)
        ax3 = fig.add_subplot(gs22[0])
        dlist = [pd.DataFrame(a.loc[x][t].agg(['max', 'min', 'median', 'mean', 'std', 'count']).round(2)).rename(columns={t: (y + " " + b)}) for x, y in ranges for a, b in iter(datah)]
        ht = dlist[0].join(dlist[1:])
        ax3.axis("off")
        ht.index = ['Maximun', 'Minimum', 'Median', 'Mean', 'std', 'Count']
        tablevalues = ht.values.tolist()
        tablevalues = [[f"{x:.2E}" if x >= 10000.0 else str(np.around(x, 2)) for x in i] for i in tablevalues]
        t1 = ax3.table(cellText=tablevalues, colLabels=(ht.columns), rowLabels=(ht.index), cellLoc="center",
          bbox=[0.0, 0, 1, 1],
          colColours=cc,
          rowColours=['#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6'])
        for (row, col), cell in t1.get_celld().items():
            if row == 0 or col == -1:
                cell.set_text_props(fontproperties=FontProperties(weight="bold"))

        for num, co in enumerate(data):
            ax1.scatter(x=(co[0]["Retention time"].values), y=(co[0][t].values), c=(dmap[num]), label=(co[1]), s=0.1, alpha=1, zorder=2, rasterized=True)

        l6 = ax11.plot((mzmla["Retention time"].values), (mzmla[TBP[l]].values), c="K", label=("MS" + str(l) + " " + TBP[l]), alpha=1, linewidth=0.2, zorder=2)
        for i in quantl:
            ax1.axvline((i[0]), 0, 1, linestyle=(i[1]), c=(i[2]), alpha=1, label=(i[3]), zorder=1)
            ax11.axvline((i[0]), 0, 1, linestyle=(i[1]), c=(i[2]), alpha=1, label=(i[3]), zorder=1)

        for i in rangesl:
            ax12.axvspan((i[0][0]), (i[0][1]), ymin=(i[1]), ymax=(i[2]), facecolor=(i[3]), edgecolor="grey")
            ax12.text(((i[0][0] + i[0][1]) / 2), (i[4]), (i[5]), horizontalalignment="center", verticalalignment="center")

        arrlistmax = []
        for b, a in enumerate(ranges):
            arrlist = [
             "", ""]
            for num, co in enumerate(datah):
                co0 = co[0].loc[a[0]]
                co1 = co[1]
                co1 = a[1] + " " + co1
                arrlist[num] = r[b].hist((co0[t].values), range=(nmin, nmax), color=(c[num][b]), alpha=1, edgecolor=(ec[b]), bins=bins, linewidth=0.5, label=co1)
                r[b].axvline(nmin, 0, 0, c=(dmap[num]), alpha=1, zorder=1, label=(legendlab[num]))

            r[b].set_xlim(nmin, nmax)
            ilabel = [(i + x) / 2 for i, x in zip(arrlist[0][1], arrlist[0][1][1:])]
            arrlistmax.append(arrlist[0][0].max())
            arrmax = max(arrlistmax)
            r[b].set_xticks(arrlist[0][1])
            if l == 2:
                if p == 1:
                    percen = zip(ilabel, arrlist[0][0] / np.sum(arrlist[0][0]) * 100, arrlist[1][0] / np.sum(arrlist[1][0]) * 100)
                else:
                    percen = zip(ilabel, arrlist[0][0] / np.sum(arrlist[0][0]) * 100)
                for i, v in enumerate(percen):
                    if v[1] != 0:
                        if b == 0:
                            heigth = 0.1
                        else:
                            heigth = 0.05
                        if len(v) == 3:
                            r[b].text((v[0]), (int(arrlist[0][0][i]) + arrmax * heigth), (str(round(v[1], 1))), color=(lgcolor[0]),
                              horizontalalignment="center",
                              verticalalignment="bottom",
                              size=8,
                              weight="bold",
                              rotation=0)
                            r[b].text((v[0]), (int(arrlist[0][0][i]) + arrmax * 0.01), (str(round(v[2], 1))), color=(lgcolor[1]),
                              horizontalalignment="center",
                              verticalalignment="bottom",
                              size=8,
                              weight="bold",
                              rotation=0)
                        else:
                            r[b].text((v[0]), (int(arrlist[0][0][i]) + arrmax * 0.01), (str(round(v[1], 1))), color=(lgcolor[0]),
                              horizontalalignment="center",
                              verticalalignment="bottom",
                              size=10,
                              weight="bold")

        ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 4), useMathText=True)
        ax11.ticklabel_format(axis="y", style="sci", scilimits=(0, 4), useMathText=True)
        Labeloffset(ax1, label=t, axis="y")
        Labeloffset(ax11, label=(l6[0].get_label()), axis="y")
        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((-3, 3))
        yfmt.set_useMathText(True)
        ax25.xaxis.set_major_formatter(yfmt)
        Labeloffset(ax25, label=t, axis="x")
        ymin, ymax = ax21.get_ylim()
        ax21.set_ylim(0, ymax * 1.2)
        ymin, ymax = ax22.get_ylim()
        ax22.set_ylim(0, ymax * 1.2)
        ax13.set_title((t + " (MS" + str(l) + ")"), fontsize="20", y=1.1)
        ax11.set_xlabel("RT", fontweight="bold")
        han, leg = ax1.get_legend_handles_labels()
        han = ax13.legend(han, leg, ncol=(len(han)), loc="center")
        if p == 1:
            han.legendHandles[5].set_sizes([11])
            han.legendHandles[6].set_sizes([11])
        else:
            han.legendHandles[5].set_sizes([11])
    for num, i in enumerate(r):
        i.set_ylabel("Counts", fontsize="10", fontweight="bold")
        i.yaxis.set_major_locator(MaxNLocator(prune="lower"))
        i.legend(fontsize=8, ncol=2)

    pp.savefig(dpi=300)
    plt.close(fig)


def bar_report(mzml, l, p, t, pp):
    id_range, aq_range, bq_range, cq_range, dq_range, rangesl, quantl = rangescal(mzml)
    maxrt = int(mzml["Retention time"].values.max())
    mzmla = mzml.loc[(mzml["ms level"] == l, ["Scan Number", "Retention time", "BPC", "TIC", t, "Scan"])]
    if l == 1:
        if p == 0:
            mzmlo = mzmla.loc[id_range]
            data = list(zip([mzmla], ["Scans"]))
            datah = list(zip([mzmlo], ["Scans"]))
            cc = ['#ffffcc', '#c2ccff', '#ffa756', '#b1d27b', '#d9544d']
            dmap = ["#0000FF"]
        if l == 2 and p == 1:
            mzml1 = mzmla[mzmla["Scan"].isna()]
            mzml2 = mzmla.dropna(0, subset=["Scan"])
            data = list(zip([mzml1, mzml2], ["Scans", "PSMs"]))
            datah = list(zip([mzmla, mzml2], ["Scans", "PSMs"]))
            cc = ['#ffffcc', '#cccc99', '#c2ccff', '#8f99cc', '#ffa756', '#cc7423',
             '#b1d27b', '#7e9f48', '#d9544d', '#a6211a']
            dmap = ["#0000FF", "#FF0000"]
    elif l == 2:
        pass
    if p == 0:
        mzmla = mzmla.dropna(0, subset=["Scan"])
        data = list(zip([mzmla], ["PSMs"]))
        datah = list(zip([mzmla], ["PSMs"]))
        cc = ['#ffffcc', '#c2ccff', '#ffa756', '#b1d27b', '#d9544d']
        dmap = ["#FF0000"]
    else:
        c = [
         [
          '#ffffcc', '#c2ccff', '#ffa756', '#b1d27b', '#d9544d'], ['#cccc99', '#8f99cc', '#cc7423', '#7e9f48', '#a6211a']]
        ec = ['#80804d', '#434d80', '#802800', '#325300', '#5a0000']
        TBP = ["", "BPC", "TIC"]
        ranges = list(zip([id_range, aq_range, bq_range, cq_range, dq_range], ['IDR', 'R1', 'R2', 'R3', 'R4']))
        legendlab = ["Scans Relative Frequency(%)", "PSMs Relative Frequency(%)"]
        lgcolor = ["#0000FF", "#FF0000"]
        fig = plt.figure(figsize=(16, 9))
        gs1 = gridspec.GridSpec(12, 1)
        gs1.update(left=0.05, right=0.52, hspace=0.0, top=0.95, bottom=0.4)
        ax1 = fig.add_subplot(gs1[2:8])
        ax11 = fig.add_subplot((gs1[8:12]), sharex=ax1)
        ax12 = fig.add_subplot((gs1[1:2]), sharex=ax1)
        ax13 = fig.add_subplot((gs1[0:1]), sharex=ax1)
        gs2 = gridspec.GridSpec(5, 1)
        gs2.update(left=0.58, right=0.98, top=0.98, bottom=0.08, hspace=0)
        ax21 = fig.add_subplot(gs2[0])
        ax21.xaxis.set_visible(False)
        ax22 = fig.add_subplot((gs2[1]), sharex=ax21)
        ax22.xaxis.set_visible(False)
        ax23 = fig.add_subplot((gs2[2]), sharey=ax22, sharex=ax22)
        ax23.xaxis.set_visible(False)
        ax24 = fig.add_subplot((gs2[3]), sharey=ax22, sharex=ax22)
        ax24.xaxis.set_visible(False)
        ax25 = fig.add_subplot((gs2[4]), sharey=ax22, sharex=ax22)
        ax1.set_xlim(0, maxrt)
        ax1.xaxis.set_visible(False)
        ax12.xaxis.set_visible(False)
        ax12.yaxis.set_visible(False)
        ax13.xaxis.set_visible(False)
        ax13.yaxis.set_visible(False)
        r = [ax21, ax22, ax23, ax24, ax25]
        gs22 = gridspec.GridSpec(1, 1)
        gs22.update(left=0.06, right=0.52, top=0.32, bottom=0.05)
        ax3 = fig.add_subplot(gs22[0])
        dlist = [pd.DataFrame(a.loc[x][t].agg(['max', 'min', 'median', 'mean', 'std', 'count']).round(2)).rename(columns={t: (y + " " + b)}) for x, y in ranges for a, b in iter(datah)]
        ht = dlist[0].join(dlist[1:])
        ax3.axis("off")
        ht.index = ['Maximun', 'Minimum', 'Median', 'Mean', 'std', 'Count']
        tablevalues = ht.values.tolist()
        tablevalues = [[f"{x:.2E}" if x >= 10000.0 else str(np.around(x, 2)) for x in i] for i in tablevalues]
        t1 = ax3.table(cellText=tablevalues, colLabels=(ht.columns), rowLabels=(ht.index), cellLoc="center",
          bbox=[0.0, 0, 1, 1],
          colColours=cc,
          rowColours=['#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6'])
        for (row, col), cell in t1.get_celld().items():
            if row == 0 or col == -1:
                cell.set_text_props(fontproperties=FontProperties(weight="bold"))

        for num, co in enumerate(data):
            ax1.scatter(x=(co[0]["Retention time"].values), y=(co[0][t].values), c=(dmap[num]), label=(co[1]), s=0.1, alpha=1, zorder=2, rasterized=True)

        l6 = ax11.plot((mzmla["Retention time"].values), (mzmla[TBP[l]].values), c="K", label=("MS" + str(l) + " " + TBP[l]), alpha=1, linewidth=0.2, zorder=2)
        ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 4))
        ax1.yaxis.set_offset_position("left")
        ax11.yaxis.set_offset_position("left")
        for i in quantl:
            ax1.axvline((i[0]), 0, 1, linestyle=(i[1]), c=(i[2]), alpha=1, label=(i[3]), zorder=1)
            ax11.axvline((i[0]), 0, 1, linestyle=(i[1]), c=(i[2]), alpha=1, label=(i[3]), zorder=1)

        for i in rangesl:
            ax12.axvspan((i[0][0]), (i[0][1]), ymin=(i[1]), ymax=(i[2]), facecolor=(i[3]), edgecolor="grey")
            ax12.text(((i[0][0] + i[0][1]) / 2), (i[4]), (i[5]), horizontalalignment="center", verticalalignment="center")

        lxaxis = []
        percenmax = []
        for b, a in enumerate(ranges):
            arrlist = [
             "", ""]
            for num, co in enumerate(datah):
                co0 = co[0].loc[a[0]]
                co1 = co[1]
                co1 = a[1] + " " + co1
                arrlist[num] = co0.groupby(t)["Scan Number"].agg("count").reset_index()
                r[b].bar((arrlist[num][t]), (arrlist[num]["Scan Number"]), color=(c[num][b]), alpha=1, edgecolor=(ec[b]), label=co1)
                r[b].axvline((min(arrlist[0][t])), 0, 0, c=(dmap[num]), alpha=1, zorder=1, label=(legendlab[num]))

            lxaxis += [min(arrlist[0][t]), max(arrlist[0][t])]
            if l == 2:
                if p == 1:
                    percen = arrlist[0].merge((arrlist[1]), how="left", left_on=t, right_on=t).fillna(value=0)
                else:
                    percen = arrlist[0].fillna(value=0)
                percenmax.append(percen.values.max())
                if b == 0:
                    arrmax = max(percenmax)
                else:
                    arrmax = max(percenmax[1:])
                for i, v in enumerate(zip(percen[t], *[percen[col] / percen[col].sum() * 100 for col in percen.columns[1:]])):
                    if v[1] != 0:
                        if len(v) == 3:
                            r[b].text((v[0]), (int(percen.iloc[(i, 1)]) + arrmax * 0.1), (str(round(v[1], 1))), color=(lgcolor[0]),
                              horizontalalignment="center",
                              verticalalignment="bottom",
                              size=8,
                              weight="bold",
                              rotation=0)
                            r[b].text((v[0]), (int(percen.iloc[(i, 1)]) + arrmax * 0.01), (str(round(v[2], 1))), color=(lgcolor[1]),
                              horizontalalignment="center",
                              verticalalignment="bottom",
                              size=8,
                              weight="bold",
                              rotation=0)
                        else:
                            r[b].text((v[0]), (int(percen.iloc[(i, 1)]) + arrmax * 0.01), (str(round(v[1], 1))), color=(lgcolor[0]),
                              horizontalalignment="center",
                              verticalalignment="bottom",
                              size=8,
                              weight="bold")

        ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 4), useMathText=True)
        ax11.ticklabel_format(axis="y", style="sci", scilimits=(0, 4), useMathText=True)
        Labeloffset(ax1, label=t, axis="y")
        Labeloffset(ax11, label=(l6[0].get_label()), axis="y")
        Labeloffset(ax25, label=t, axis="x")
        ymin, ymax = ax21.get_ylim()
        ax21.set_ylim(0, ymax * 1.2)
        ymin2, ymax2 = ax22.get_ylim()
        ax22.set_ylim(0, ymax2 * 1.2)
        if max(lxaxis) > 35:
            lxaxis = list(range(int(min(lxaxis)), int(max(lxaxis)) + 1, 2))
        else:
            lxaxis = list(range(int(min(lxaxis)), int(max(lxaxis)) + 1))
        ax24.set_xticks(lxaxis)
        ax13.set_title((t + " (MS" + str(l) + ")"), fontsize="20", y=1.1)
        ax11.set_xlabel("RT", fontweight="bold")
        han, leg = ax1.get_legend_handles_labels()
        han = ax13.legend(han, leg, ncol=(len(han)), loc="center")
        if p == 1:
            han.legendHandles[5].set_sizes([11])
            han.legendHandles[6].set_sizes([11])
        else:
            han.legendHandles[5].set_sizes([11])
    for num, i in enumerate(r):
        i.set_ylabel("Frequency", fontsize="10", fontweight="bold")
        i.yaxis.set_major_locator(MaxNLocator(prune="lower"))
        i.legend(fontsize=8, ncol=2)

    pp.savefig(dpi=300)
    plt.close(fig)


def main_report(mzml, pp, label):
    id_range, aq_range, bq_range, cq_range, dq_range, rangesl, quantl = rangescal(mzml)
    maxrt = int(mzml["Retention time"].values.max())
    cc = ['#ffffcc', '#c2ccff', '#ffa756', '#b1d27b', '#d9544d']
    TBP = list(zip(["BPC", "TIC"], [1, 2]))
    ranges = list(zip([id_range, aq_range, bq_range, cq_range, dq_range], ['IDR', 'R1', 'R2', 'R3', 'R4']))
    fig = plt.figure(figsize=(16, 9))
    gs1 = gridspec.GridSpec(12, 1)
    gs1.update(left=0.05, right=0.52, hspace=0.0, top=0.95, bottom=0.4)
    ax1 = fig.add_subplot(gs1[2:7])
    ax11 = fig.add_subplot((gs1[7:12]), sharex=ax1)
    ax12 = fig.add_subplot((gs1[1:2]), sharex=ax1)
    ax13 = fig.add_subplot((gs1[0:1]), sharex=ax1)
    gs2 = gridspec.GridSpec(1, 1)
    gs2.update(left=0.65, right=0.98, top=0.98, bottom=0.02, hspace=0.28)
    ax2 = fig.add_subplot(gs2[0])
    ax1.set_xlim(0, maxrt)
    ax1.xaxis.set_visible(False)
    ax12.xaxis.set_visible(False)
    ax12.yaxis.set_visible(False)
    ax13.xaxis.set_visible(False)
    ax13.yaxis.set_visible(False)
    lir = [ax1, ax11]
    gs22 = gridspec.GridSpec(1, 1)
    gs22.update(left=0.1, right=0.52, top=0.34, bottom=0.02)
    ax3 = fig.add_subplot(gs22[0])
    llist = [
     "", ""]
    for b, i in enumerate(TBP):
        mzml1 = mzml.loc[mzml["ms level"] == i[1]]
        llist[b] = lir[b].plot((mzml1["Retention time"].values), (mzml1[i[0]].values), c="K", label=("MS" + str(b + 1) + " " + i[0]), alpha=1, linewidth=0.2, zorder=2)

    for i in quantl:
        ax1.axvline((i[0]), 0, 1, linestyle=(i[1]), c=(i[2]), alpha=1, label=(i[3]), zorder=1)
        ax11.axvline((i[0]), 0, 1, linestyle=(i[1]), c=(i[2]), alpha=1, label=(i[3]), zorder=1)

    for i in rangesl:
        ax12.axvspan((i[0][0]), (i[0][1]), ymin=(i[1]), ymax=(i[2]), facecolor=(i[3]), edgecolor="grey")
        ax12.text(((i[0][0] + i[0][1]) / 2), (i[4]), (i[5]), horizontalalignment="center", verticalalignment="center", fontweight="bold")

    table1 = []
    for i in ranges:
        mzml1 = mzml.loc[i[0]]
        table11 = list(zip([len(mzml1[mzml1["ms level"] == 1])], [
         len(mzml1[mzml1["ms level"] == 2])], [
         len(mzml1.dropna(subset=["Scan"]))], [
         len(mzml1["Modified Sequence"].unique())], [
         len(mzml1.Protein.unique())], [
         len(mzml1[(mzml1["Intensity jumps"] > 10) | (mzml1["Intensity jumps"] < 0.1)])], [
         str(round(time_atr(mzml1)[0], 1))], [
         str(round(time_atr(mzml1)[1], 1))], [
         str(round(time_atr(mzml1)[2], 1))]))
        table1 += table11

    ms1 = ['Injection time', 'Intensity', 'm/z count', 'TopN', 'ID efficiency [%]',
     'Median Intensity']
    if label != 4:
        ms2 = ['Injection time', 'Intensity', 'm/z count', 'Precursor m/z', 'Hyperscore',
         'Delta Mass [ppm]', 'Matched ions [%]', 'Charge', 'Missed cleavages',
         'Redundancy', 'Median Intensity', 'Reporter count', 'Median reporter int']
        ind = ['Injection time (MS1)', 'Intensity (MS1)', 'm/z count (MS1)', 'TopN',
         'ID efficiency [%]', 'Median Intensity (MS1)', 'Injection time (MS2)',
         'Intensity (MS2)', 'm/z count (MS2)', 'Precursor m/z', 'Hyperscore',
         'Delta Mass [ppm]', 'Matched ions [%]', 'Charge', 'Missed cleavages',
         'Redundancy', 'Median Intensity (MS2)', 'Reporter count', 'Median reporter int']
        loc1 = ['Redundancy', 'Intensity (MS1)', 'Median Intensity (MS1)', 'Injection time (MS1)',
         'm/z count (MS1)', 'Intensity (MS2)', 'Median Intensity (MS2)',
         'Injection time (MS2)', 'm/z count (MS2)', 'TopN', 'Charge',
         'Precursor m/z', 'Reporter count', 'Median reporter int', 'ID efficiency [%]',
         'Delta Mass [ppm]', 'Hyperscore', 'Matched ions [%]', 'Missed cleavages']
    else:
        ms2 = [
         'Injection time', 'Intensity', 'm/z count', 'Precursor m/z',
         'Hyperscore', 'Delta Mass [ppm]', 'Matched ions [%]', 'Charge',
         'Missed cleavages', 'Redundancy', 'Median Intensity']
        ind = ['Injection time (MS1)', 'Intensity (MS1)', 'm/z count (MS1)', 'TopN',
         'ID efficiency [%]', 'Median Intensity (MS1)', 'Injection time (MS2)',
         'Intensity (MS2)', 'm/z count (MS2)', 'Precursor m/z', 'Hyperscore',
         'Delta Mass [ppm]', 'Matched ions [%]', 'Charge', 'Missed cleavages',
         'Redundancy', 'Median Intensity (MS2)']
        loc1 = ['Redundancy', 'Intensity (MS1)', 'Median Intensity (MS1)', 'Injection time (MS1)',
         'm/z count (MS1)', 'Intensity (MS2)', 'Median Intensity (MS2)',
         'Injection time (MS2)', 'm/z count (MS2)', 'TopN', 'Charge',
         'Precursor m/z', 'ID efficiency [%]', 'Delta Mass [ppm]', 'Hyperscore',
         'Matched ions [%]', 'Missed cleavages']
    table1 = pd.DataFrame(table1).T
    table1.columns = ['IDR', 'R1', 'R2', 'R3', 'R4']
    table1.index = ['MS1_Scans', 'MS2_Scans', 'PSMs', 'Peptides', 'Proteins', 'Intensity jumps',
     'RT(min)', 'RT(max)', 'RT(diff)']
    t1 = ax3.table(cellText=(table1.values), colLabels=(table1.columns), rowLabels=(table1.index), cellLoc="center",
      bbox=[0.0, 0, 1, 1],
      colColours=cc,
      rowColours=['#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6',
     '#d8dcd6', '#d8dcd6'])
    ax3.axis("off")
    mzml = mzml.loc[id_range]
    table2 = mzml.loc[(mzml["ms level"] == 1, ms1)].agg(["max", "min", "median", "mean"]).T
    table2 = table2.append(mzml.loc[(mzml["ms level"] == 2, ms2)].agg(["max", "min", "median", "mean"]).T)
    table2.index = ind
    table2 = table2.loc[loc1, :]
    table2.columns = ["IDR Maximun", "IDR Minimum", "IDR Median", "IDR Mean"]
    table2values = table2.values.tolist()
    table2values = [[f"{x:.2E}" if x >= 10000.0 else str(np.around(x, 2)) for x in i] for i in table2values]
    rowcolours = ["#d8dcd6" for i in table2.index]
    t2 = ax2.table(cellText=table2values, colLabels=(table2.columns), rowLabels=(table2.index), cellLoc="center",
      bbox=[0.0, 0, 1, 1],
      rowColours=rowcolours,
      colColours=["#d8dcd6", "#d8dcd6", "#d8dcd6", "#d8dcd6"])
    ax2.axis("off")
    tl = [t1, t2]
    for i in tl:
        for (row, col), cell in i.get_celld().items():
            if row == 0 or col == -1:
                cell.set_text_props(fontproperties=FontProperties(weight="bold"))

    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 4), useMathText=True)
    ax11.ticklabel_format(axis="y", style="sci", scilimits=(0, 4), useMathText=True)
    Labeloffset(ax1, label=(llist[0][0].get_label()), axis="y")
    Labeloffset(ax11, label=(llist[1][0].get_label()), axis="y")
    ax13.set_title("Summary", fontsize="20", y=1.1)
    ax11.set_xlabel("RT", fontweight="bold")
    han, leg = ax1.get_legend_handles_labels()
    han = han[1:]
    leg = leg[1:]
    ax13.legend(han, leg, ncol=(len(han)), loc="center")
    pp.savefig()
    plt.close(fig)


def sampling_report(mzml, pp):
    id_range, aq_range, bq_range, cq_range, dq_range, rangesl, quantl = rangescal(mzml)
    maxrt = int(mzml["Retention time"].values.max())
    TBP = ["Acummulated redundancy", "Redundancy", "Resolution range"]
    ranges = list(zip([id_range, aq_range, bq_range, cq_range, dq_range], ['IDR', 'R1', 'R2', 'R3', 'R4']))
    mzml2 = mzml[mzml["ms level"] == 2]
    mzml1 = mzml.dropna(0, subset=["Scan"])
    fig = plt.figure(figsize=(16, 9))
    gs1 = gridspec.GridSpec(20, 1)
    gs1.update(left=0.05, right=0.52, hspace=0.0, top=0.95, bottom=0.08)
    ax1 = fig.add_subplot(gs1[2:7])
    ax11 = fig.add_subplot((gs1[7:12]), sharex=ax1)
    ax14 = fig.add_subplot((gs1[12:17]), sharex=ax1)
    ax15 = fig.add_subplot((gs1[17:20]), sharex=ax1)
    ax12 = fig.add_subplot((gs1[1:2]), sharex=ax1)
    ax13 = fig.add_subplot((gs1[0:1]), sharex=ax1)
    gs2 = gridspec.GridSpec(1, 1)
    gs2.update(left=0.65, right=0.98, top=0.98, bottom=0.08, hspace=0.28)
    ax2 = fig.add_subplot(gs2[0])
    ax1.set_xlim(0, maxrt)
    ax1.xaxis.set_visible(False)
    ax11.xaxis.set_visible(False)
    ax14.xaxis.set_visible(False)
    ax12.xaxis.set_visible(False)
    ax12.yaxis.set_visible(False)
    ax13.xaxis.set_visible(False)
    ax13.yaxis.set_visible(False)
    lir = [ax1, ax11, ax14]
    dmap = [["#FF0000", "#922428", "#922428"], ["PSMs", "Peptides", "Peptides"]]
    llist = [
     "", "", ""]
    for b, i in enumerate(TBP):
        llist[b] = lir[b].scatter(x=(mzml1["Retention time"].values), y=(mzml1[i].values), c=(dmap[0][b]), label=(dmap[1][b]), alpha=1, zorder=2, rasterized=True, s=0.1)

    ax15.plot((mzml2["Retention time"].values), (mzml2["Intensity"].values), c="K", label="MS2 TIC", alpha=1, linewidth=0.2, zorder=2)
    for i in quantl:
        ax1.axvline((i[0]), 0, 1, linestyle=(i[1]), c=(i[2]), alpha=1, label=(i[3]), zorder=1)
        ax11.axvline((i[0]), 0, 1, linestyle=(i[1]), c=(i[2]), alpha=1, label=(i[3]), zorder=1)
        ax14.axvline((i[0]), 0, 1, linestyle=(i[1]), c=(i[2]), alpha=1, label=(i[3]), zorder=1)
        ax15.axvline((i[0]), 0, 1, linestyle=(i[1]), c=(i[2]), alpha=1, label=(i[3]), zorder=1)

    for i in rangesl:
        ax12.axvspan((i[0][0]), (i[0][1]), ymin=(i[1]), ymax=(i[2]), facecolor=(i[3]), edgecolor="grey")
        ax12.text(((i[0][0] + i[0][1]) / 2), (i[4]), (i[5]), horizontalalignment="center", verticalalignment="center", fontweight="bold")

    table1 = pd.DataFrame()
    for i in ranges:
        table2 = pd.DataFrame(mzml1.loc[(i[0], ["Redundancy", "Resolution range"])].agg(['max', 'min', 'median', 'mean', 'count']).stack())
        table2.columns = [i[1]]
        table1 = pd.concat([table1, table2], axis=1)

    table1 = table1.reset_index().sort_values(by="level_1")
    table1["in"] = table1["level_1"] + [" "] + table1["level_0"]
    table1 = table1.set_index("in")
    table1 = table1[['IDR', 'R1', 'R2', 'R3', 'R4']]
    table1values = table1.values.tolist()
    table1values = [[f"{x:.2E}" if x >= 10000.0 else str(np.around(x, 2)) for x in i] for i in table1values]
    rowcolours = ['#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6', '#d8dcd6', '#a5a9a3', '#a5a9a3',
     '#a5a9a3', '#a5a9a3', '#a5a9a3']
    t1 = ax2.table(cellText=table1values, colLabels=(table1.columns), rowLabels=(table1.index), cellLoc="center",
      bbox=[0.0, 0, 1, 1],
      rowColours=rowcolours,
      colColours=['#ffffcc', '#c2ccff', '#ffa756', '#b1d27b', '#d9544d'])
    ax2.axis("off")
    for (row, col), cell in t1.get_celld().items():
        if row == 0 or col == -1:
            cell.set_text_props(fontproperties=FontProperties(weight="bold"))

    ax13.set_title("Redundancy & Resolution range (Part 1)", fontsize="20", y=1.1)
    ax11.set_ylabel("Redundancy", fontweight="bold")
    ax14.set_ylabel("Resolution range", fontweight="bold")
    ax15.ticklabel_format(axis="y", style="sci", scilimits=(0, 4), useMathText=True)
    Labeloffset(ax15, label="MS2 TIC", axis="y")
    ax15.set_xlabel("RT", fontweight="bold")
    ax1.set_ylabel("Acummulated redundancy", fontweight="bold")
    han, leg = ax1.get_legend_handles_labels()
    han1, leg1 = ax11.get_legend_handles_labels()
    han.append(han1[-1])
    leg.append(leg1[-1])
    han = ax13.legend(han, leg, ncol=(len(han)), loc="center")
    han.legendHandles[5].set_sizes([11])
    han.legendHandles[6].set_sizes([11])
    pp.savefig(dpi=300)
    plt.close(fig)
    fig = plt.figure(figsize=(16, 9))
    gs1 = gridspec.GridSpec(5, 1)
    gs1.update(left=0.05, right=0.48, top=0.95, bottom=0.08, hspace=0)
    ax11 = fig.add_subplot(gs1[0])
    ax11.xaxis.set_visible(False)
    ax12 = fig.add_subplot((gs1[1]), sharex=ax11)
    ax12.xaxis.set_visible(False)
    ax13 = fig.add_subplot((gs1[2]), sharey=ax12, sharex=ax11)
    ax13.xaxis.set_visible(False)
    ax14 = fig.add_subplot((gs1[3]), sharey=ax12, sharex=ax11)
    ax14.xaxis.set_visible(False)
    ax15 = fig.add_subplot((gs1[4]), sharey=ax12, sharex=ax11)
    gs2 = gridspec.GridSpec(5, 1)
    gs2.update(left=0.55, right=0.98, top=0.95, bottom=0.08, hspace=0)
    ax21 = fig.add_subplot(gs2[0])
    ax21.xaxis.set_visible(False)
    ax22 = fig.add_subplot((gs2[1]), sharex=ax21)
    ax22.xaxis.set_visible(False)
    ax23 = fig.add_subplot((gs2[2]), sharey=ax22, sharex=ax21)
    ax23.xaxis.set_visible(False)
    ax24 = fig.add_subplot((gs2[3]), sharey=ax22, sharex=ax21)
    ax24.xaxis.set_visible(False)
    ax25 = fig.add_subplot((gs2[4]), sharey=ax22, sharex=ax21)
    r1 = [
     [
      ax11, ax12, ax13, ax14, ax15], [ax21, ax22, ax23, ax24, ax25]]
    dat = [mzml["Redundancy"], mzml["Resolution range"]]
    c = [
     '#ffffcc', '#c2ccff', '#ffa756', '#b1d27b', '#d9544d']
    ec = ['#80804d', '#434d80', '#802800', '#325300', '#5a0000']
    for x, y in enumerate(r1):
        nmax = round(outQ(dat[x].dropna()).max(), 1)
        nmin = round(outQ(dat[x].dropna()).min(), 1)
        bins, nmax, nmin = bins_range(nmax, nmin)
        for b, a in enumerate(ranges):
            co = a[1] + " " + "Peptides"
            arrlist = y[b].hist((dat[x][a[0]].values), range=(nmin, nmax), color=(c[b]), alpha=1, edgecolor=(ec[b]), bins=bins, linewidth=0.5, label=co)
            y[b].set_xlim(nmin, nmax)
            y[b].set_xticks(arrlist[1])
            yfmt = ScalarFormatterForceFormat()
            yfmt.set_powerlimits((-3, 3))
            y[b].xaxis.set_major_formatter(yfmt)
            y[b].axvline(nmin, 0, 0, c="#922428", alpha=1, zorder=1, label="Peptides Relative Frequency(%)")
            ilabel = [(i + x) / 2 for i, x in zip(arrlist[1], arrlist[1][1:])]
            percen = zip(ilabel, arrlist[0] / np.sum(arrlist[0]) * 100)
            for i, v in enumerate(percen):
                y[b].text((v[0]), (int(arrlist[0][i]) + arrlist[0].max() * 0.01), (str(round(v[1], 1))), color="#922428",
                  horizontalalignment="center",
                  verticalalignment="bottom",
                  size=10,
                  weight="bold")

    ax11.set_title("Redundancy (Part 2)", fontsize="20", y=1.05)
    ax21.set_title("Resolution range (Part 2)", fontsize="20", y=1.05)
    ax15.set_xlabel("Redundancy", fontsize="10", fontweight="bold")
    ax25.set_xlabel("Resolution range", fontsize="10", fontweight="bold")
    ymin, ymax = ax11.get_ylim()
    ax11.set_ylim(0, ymax * 1.2)
    ymin2, ymax2 = ax12.get_ylim()
    ax12.set_ylim(0, ymax2 * 1.2)
    ymin1, ymax1 = ax21.get_ylim()
    ax21.set_ylim(0, ymax1 * 1.2)
    ymin12, ymax12 = ax22.get_ylim()
    ax22.set_ylim(0, ymax12 * 1.2)
    for x in r1:
        for num, i in enumerate(x):
            i.set_ylabel("Counts", fontsize="10", fontweight="bold")
            i.yaxis.set_major_locator(MaxNLocator(prune="lower"))
            i.legend(fontsize=8, ncol=2)

    pp.savefig()
    plt.close(fig)