####################################################################
# --- This script takes the plotter config as its input. --------- #
# --- It took the last bin content from 'evtCutFlow' histogram --- #
# --- to estimate the bkg contribution in the SR along with -------#   
# --- a matrix of signal yields for (6x6) coupling points ---------#
####################################################################

import os
import sys
import copy
import yaml
import ROOT
import math
from collections import defaultdict
import logging
import argparse
from prettytable import PrettyTable
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def getHeatmap(data, row_labels, col_labels, ax=None,
               cbar_kw={}, cbarlabel="", **kwargs):
    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)
    ax.set_xlabel(r"$Y_{l}$")
    ax.set_ylabel(r"$Y_{q}$")

    # Let the horizontal axes labeling appear on top.
    #ax.tick_params(top=True, bottom=False,
    #               labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    #ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

def computeSignificance(nSig, nBkg, isLog=False):
    signf = math.sqrt(2*(nSig+nBkg)*math.log(1+(nSig/nBkg)) - 2*nSig) if isLog else nSig/math.sqrt(nSig+nBkg)
    return signf

def main():
    parser = argparse.ArgumentParser(description='JobCardMaker')
    parser.add_argument('--Process', type=str, action='store', required=True, help='XZ or HZ')
    parser.add_argument('--MX', type=str, action='store', required=True, help='')
    parser.add_argument('--MH', type=str, action='store', required=True, help='')
    parser.add_argument('--yaml', action='store', required=True, type=str, help='Name of the config')

    args = parser.parse_args()
    Yl = [0.001,0.003,0.005,0.007,0.009,0.01]
    Yq = Yl
    tag = args.Process+'_'+str(args.MX)+'_'+str(args.MH)
    filetag = str(args.yaml).split('.')[0]
    couplList = []
    for yl in Yl:
        for yq in Yq:
            couplList.append((yl,yq))
    print(f'List of couplings (Yl,Yq) : {couplList}')

    with open(os.path.join('Yamls',args.yaml), 'r') as config:
        configDict = yaml.safe_load(config)

    lumi     = configDict.get('Lumi')
    InputDir = configDict.get('InDir')
    plotDir  = os.path.join(os.getcwd(),'Plots_'+filetag)

    if os.path.exists(plotDir):
        print(f'{plotDir} exists. It will be overwritten!!!')
    else:
        print(f'Creating plot dir: {plotDir}')
        os.mkdir(plotDir)

    t0 = PrettyTable(['Process', f'Yield at {lumi/1000} fb-1'])
    t0.title = f'{tag}'

    sampleDict   = configDict.get('Samples')
    plotExt      = configDict.get('plotExtension')
    handle = args.Process+'_'+str(args.MX)+'_'+str(args.MH)
    logf   = open(os.path.join(plotDir, str(args.yaml).split('.')[0]+'_Significance.log'), 'w') 
    logf.write('Yields at {lumi/1000} fb-1 \n')
    yieldlog = open(os.path.join(plotDir, str(args.yaml).split('.')[0]+'_Yield.log'), 'w')
    nBkgs  = 0
    nSignal1dArray = [0]*36
    for sampleKey, sampleValDict in sampleDict.items():
        infile = os.path.join(InputDir,sampleValDict.get('file'))
        logging.info(f'Sample Name : {sampleKey}, file : {infile}')
        xsec   = sampleValDict.get('cross-section')
        label  = sampleValDict.get('label')
        rfile  = ROOT.TFile(infile,'READ')
        if rfile.IsZombie():
            logger.warning(f'W A R N I N G ::: {infile} is a Zombie !!!')
            logf.write(f'can not open {infile}')
            continue
        evtCutFlow_hist  = copy.deepcopy(rfile.Get('evtCutFlow'))
        if (evtCutFlow_hist.Integral() > 0):
            evtCutFlow_hist.Scale(lumi*xsec/evtCutFlow_hist.GetBinContent(evtCutFlow_hist.FindBin(0)))
        else:
            print('W A R N I N G ::: evtCutFlow histogram cant be scaled by lumi as integral = 0')

        endBin = evtCutFlow_hist.FindBin(evtCutFlow_hist.GetNbinsX()) - 1
        endBinContent = evtCutFlow_hist.GetBinContent(endBin)

        logf.write(f'{infile} : --->  {endBinContent} \n')
        t0.add_row([sampleKey, endBinContent])
        for i, coup in enumerate(couplList):
            key = tag+'_Yl_'+str(coup[0])+'_Yq_'+str(coup[1])
            #print(f'{sampleKey}:{key}')
            if sampleKey == key:
                nSignal1dArray[i] = endBinContent
            #else:
            #    print(f'{coup} couplings not found in yaml')
        
        if not 'XZ' in sampleKey:
            #print(f'{sampleKey} : {endBinContent}')
            #logf.write(f'{sampleKey} : {endBinContent} \n')
            nBkgs += endBinContent
    
    t0.add_row(['Total Bkg', nBkgs])
    yieldlog.write(str(t0)+'\n')
    significance1d    = np.array([computeSignificance(nSig, nBkgs, isLog=False) for nSig in nSignal1dArray])
    significance1dLog = np.array([computeSignificance(nSig, nBkgs, isLog=True) for nSig in nSignal1dArray])
    print(significance1d.shape)
    significance1d.resize(6,6)
    print(significance1d.shape)
    significance1dLog.resize(6,6)

    significance1dSyst = np.array([computeSignificance(nSig, nBkgs + 0.01*nBkgs*nBkgs, isLog=False) for nSig in nSignal1dArray])
    significance1dSyst.resize(6,6)

    significanceMatrix     = np.flipud(significance1d)
    significanceMatrixLog  = np.flipud(significance1dLog)
    significanceMatrixSyst = np.flipud(significance1dSyst)

    np.array(nSignal1dArray).resize(6,6)
    nSig2d = np.flipud(nSignal1dArray)

    logf.write('Signal yields at {lumi/1000} fb-1 \n')
    logf.write(str(nSig2d)+'\n')

    print(f'by S/Sqrt(S+B) : \n {significanceMatrix}')
    print(f'by log formula : \n {significanceMatrixLog}')
    print(f'with systematics: \n {significanceMatrixSyst}')
    logf.write(f'S/Sqrt(S+B) Significance of {filetag} Analysis \n')
    logf.write(str(significanceMatrix)+'\n')
    logf.write(f'Log Significance of {filetag} Analysis \n')
    logf.write(str(significanceMatrixLog)+'\n')
    logf.write(f'Significance of {filetag} Analysis with 10% systematics\n')
    logf.write(str(significanceMatrixSyst)+'\n')

    fig, ax = plt.subplots()
    im = ax.imshow(significanceMatrix)
    im, cbar = getHeatmap(significanceMatrix, Yq[::-1], Yl, ax=ax, cmap="YlGn", cbarlabel="Scale")
    texts = annotate_heatmap(im, valfmt="{x:.2f}")
    fig.tight_layout()
    plt.savefig(os.path.join(plotDir,handle+'_NotLog_2dSignificance'+'.'+plotExt), format=plotExt)
    plt.clf()
    
    fig, ax = plt.subplots()
    im = ax.imshow(significanceMatrixLog)
    im, cbar = getHeatmap(significanceMatrixLog, Yq[::-1], Yl, ax=ax, cmap="YlGn", cbarlabel="Scale")
    texts = annotate_heatmap(im, valfmt="{x:.2f}")
    fig.tight_layout()
    plt.savefig(os.path.join(plotDir,handle+'_Log_2dSignificance'+'.'+plotExt), format=plotExt)
    plt.clf()

    fig, ax = plt.subplots()
    im = ax.imshow(significanceMatrixSyst)
    im, cbar = getHeatmap(significanceMatrixSyst, Yq[::-1], Yl, ax=ax, cmap="YlGn", cbarlabel="Scale")
    texts = annotate_heatmap(im, valfmt="{x:.2f}")
    fig.tight_layout()
    plt.savefig(os.path.join(plotDir,handle+'_With10pSyst_2dSignificance'+'.'+plotExt), format=plotExt)
    plt.clf()

    significanceYl0p001 = list(np.flip(significanceMatrix[:,0]))
    significanceYl0p01  = list(np.flip(significanceMatrix[:,-1]))
    significanceYl0p001Syst = list(np.flip(significanceMatrixSyst[:,0]))
    significanceYl0p01Syst  = list(np.flip(significanceMatrixSyst[:,-1]))

    print(significanceYl0p001, significanceYl0p01)
    print(significanceYl0p001Syst, significanceYl0p01Syst)

    Yq = ['0.001','0.003','0.005','0.007','0.009','0.01']
    sns.set_style("white")
    plt.figure(figsize=(11, 8.5))
    fig, ax = plt.subplots(1)
    ax.plot(Yq, significanceYl0p001, label=r"$Y_{\ell} = 0.001$", color='green', lw=2)
    ax.plot(Yq, significanceYl0p01, color='green', linestyle='--', label=r"$Y_{\ell} = 0.01$", lw=2)
    ax.plot(Yq, significanceYl0p001Syst, label=r"$Y_{\ell} = 0.001$", color='red', lw=2)
    ax.plot(Yq, significanceYl0p01Syst, color='red', linestyle='--', label=r"$Y_{\ell} = 0.01$", lw=2)

    ax.fill_between(Yq, significanceYl0p001, significanceYl0p01, interpolate=True, facecolor='green', hatch='\\', alpha=0.7, label='Nominal')
    ax.fill_between(Yq, significanceYl0p001Syst, significanceYl0p01Syst, interpolate=True, facecolor='red', hatch='/', alpha=0.4, label='With 10% Systematics')

    ax.set_xlabel(r"$Y_{q}$",labelpad=3.0,fontsize=14)
    ax.set_ylabel('Estimated significance',fontsize=14)
    #ax.set_xticks(20)                                                                                                                             
    #ax.set_yticks(20)                                                                                                                             
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    #ax.grid()
    '''
    custom_lines = [Line2D([0], [0], color='green', alpha=0.7, lw=8),
                    Line2D([0], [0], color='red', alpha=0.4, lw=8),
                    Line2D([0], [0], color='black', alpha=0.8, lw=3, linestyle='-'),
                    Line2D([0], [0], color='black', alpha=0.8, lw=3, linestyle='--')]
    leg = ax.legend(loc='best',
                    fontsize=15,
                    fancybox=True,
                    facecolor='None')
    leg.legendHandles[0].set_color('green')
    leg.legendHandles[0].set_alpha(0.4)
    leg.legendHandles[1].set_color('red')
    leg.legendHandles[1].set_alpha(0.4)
    custom_lines = [Line2D([0], [0], color='green', alpha=0.7, lw=8),
                    Line2D([0], [0], color='red', alpha=0.4, lw=8)]
    '''
    custom_lines = [mpatches.Patch(facecolor="green",alpha=0.7,hatch='\\\\',label='Nominal'),
                    mpatches.Patch(facecolor="red",alpha=0.4,hatch='//',label='10 % systematics')]

    ax.legend(handles=custom_lines,
              #title=r'$Y_{#ell} = 0.001 \to 0.01$',
              loc='best',
              fontsize=15,
              fancybox=True,
              facecolor='None')
    plt.tight_layout()
    plt.savefig(os.path.join(plotDir,handle+'_1dSignificance.png'), dpi=300)

if __name__ == "__main__":
    main()
