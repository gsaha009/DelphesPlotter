###################################################
#            An easy Plotter Script               #
#             Author : Gourab Saha                #
# Desription: It takes an yaml as input. Keep the #
# yaml files in 'Yamls' dir. This script ables to #
# print the yield tables along with the plots.    #
################################################### 

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
#from memory_profiler import profile

def getListOfHistograms(inputDir,sampleDict):
    logger = logging.getLogger("plotLog")
    for i,(key,val) in enumerate(sampleDict.items()):
        if (i > 0):
            break
        file_ = os.path.join(inputDir, val.get('file'))
        if not os.path.exists(file_):
            logger.warning(f'{file_} doesnt exist')
        infile = ROOT.TFile(file_,'READ')
        histList  = [hkey.GetName() for hkey in infile.GetListOfKeys() 
                     if not any(['ObjectSelection' in hkey.GetName(),'LepPairSelection' in hkey.GetName(), 'CutFlow' in hkey.GetName()])]
        titleList = [hkey.GetTitle() for hkey in infile.GetListOfKeys() 
                     if not any(['ObjectSelection' in hkey.GetName(),'LepPairSelection' in hkey.GetName(), 'CutFlow' in hkey.GetName()])]
    return histList, titleList

def getCanvas(canvasName, legpos, mergin, canvSize, legendTxtSize=0.04, logy=False):
    canv = ROOT.TCanvas(canvasName, "Plot", canvSize[0], canvSize[1])
    legend = ROOT.TLegend(legpos[0],
                          legpos[1],
                          legpos[2],
                          legpos[3])
    legend.SetName(canvasName+'_leg')
    legend.SetBorderSize(1)
    legend.SetTextSize(legendTxtSize)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)

    ROOT.gPad.SetLeftMargin(mergin[0])
    ROOT.gPad.SetRightMargin(mergin[1])
    ROOT.gPad.SetTopMargin(mergin[2])
    ROOT.gPad.SetBottomMargin(mergin[3])
    if logy:
        ROOT.gPad.SetLogy()
    #ROOT.gPad.SetGrid(1)
    
    ROOT.gStyle.SetTitleFontSize(0.05)
    ROOT.gStyle.SetOptStat(0)
    
    return canv, legend

def getHistPerLabel(plotlist):
    for ip, plot in enumerate(plotlist):
        if ip == 0:
            plot_copy = copy.deepcopy(plot)
        else:
            plot_copy.Add(plot)
    if (plot_copy.Integral() > 0):
        plot_copy.Scale(1/plot_copy.Integral())
    return plot_copy

def getEfficiencies(ibin, evtCutFlow_hist):
    if ibin == 0:
        abs_eff = 1.0
        rel_eff = 1.0
    else:
        den_abs = evtCutFlow_hist.GetBinContent(evtCutFlow_hist.FindBin(0))
        den_rel = evtCutFlow_hist.GetBinContent(evtCutFlow_hist.FindBin(ibin-1))
        num     = evtCutFlow_hist.GetBinContent(evtCutFlow_hist.FindBin(ibin))

        abs_eff = num/den_abs if den_abs > 0 else math.nan
        rel_eff = num/den_rel if den_rel > 0 else math.nan

    return [abs_eff, rel_eff]


def getYields(sampleDict, InputDir, lumi):
    logger = logging.getLogger("plotLog")
    tables = []
    for sampleKey, sampleValDict in sampleDict.items():
        infile = os.path.join(InputDir,sampleValDict.get('file'))
        logging.info(f'Sample Name : {sampleKey}, file : {infile}')
        if not os.path.isfile(infile):
            print(f'W A R N I N G :: -----> {infile} not found !!!')
            continue
        xsec   = sampleValDict.get('cross-section')
        label  = sampleValDict.get('label')
        rfile  = ROOT.TFile(infile,'READ')
        if rfile.IsZombie():
            logger.warning(f'{infile} is a Zombie !!!')
            continue
            
        t0 = PrettyTable(['Process', 'Yield : Unweighted', f'Yield : Lumi-weighted', 'Relative Efficiency', 'Absolute Efficiency','effective xsec (fb)'])
        t0.title = f'Process : {sampleKey}, X-sec : {xsec*1000} fb, Lumi : {lumi/1000} fb-1'
        evtCutFlow_hist      = copy.deepcopy(rfile.Get('evtCutFlow'))
        evtCutFlow_hist_copy = copy.deepcopy(evtCutFlow_hist)
        if (evtCutFlow_hist_copy.Integral() > 0):
            evtCutFlow_hist_copy.Scale(lumi*xsec/evtCutFlow_hist_copy.GetBinContent(evtCutFlow_hist_copy.FindBin(0)))
        else:
            print('evtCutFlow histogram cant be scaled by lumi as integral = 0')
        

        for i in range(evtCutFlow_hist.GetNbinsX()):
            effList = getEfficiencies(i, evtCutFlow_hist)
            t0.add_row([evtCutFlow_hist.GetXaxis().GetBinLabel(evtCutFlow_hist.FindBin(i)), 
                        round(evtCutFlow_hist.GetBinContent(evtCutFlow_hist.FindBin(i)),2), 
                        round(evtCutFlow_hist_copy.GetBinContent(evtCutFlow_hist_copy.FindBin(i)),2), 
                        '%.3f'%effList[1],
                        '%.3f'%effList[0],
                        '%.3f'%(effList[0]*float(xsec)*1000)])
        rfile.Close()
        tables.append(t0)
    return tables

def main():
    logging.basicConfig(level   = logging.DEBUG,
                        format  = '%(asctime)s - %(levelname)s - %(message)s',
                        datefmt = '%m/%d/%Y %H:%M:%S')
    logger = logging.getLogger("plotLog")
    logger.setLevel(logging.DEBUG)    

    parser = argparse.ArgumentParser(description='Histogram plotter')    
    parser.add_argument('--yaml', action='store', required=True, type=str, help='Name of the config')
    parser.add_argument('--yieldOnly', action='store_true',required=False, default=False, help='Print the yield tables')

    args = parser.parse_args()

    pwd = os.getcwd()
    logger.info(f'Current working directory : {pwd}')

    with open(os.path.join('Yamls',args.yaml), 'r') as config:
        configDict = yaml.safe_load(config)

    lumi     = configDict.get('Lumi')
    InputDir = configDict.get('InDir') 
    plotDir  = os.path.join(pwd,'Plots_'+str(args.yaml).split('.')[0])

    if os.path.exists(plotDir):
        logger.info(f'{plotDir} exists. It will be overwritten!!!')
    else:
        logger.info(f'Creating plot dir: {plotDir}')
        os.mkdir(plotDir)

    sampleDict   = configDict.get('Samples')
    plotExt      = configDict.get('plotExtension')
    canvSize     = configDict.get('CanvSize')
    legendPos    = configDict.get('legendPos')
    mergin       = configDict.get('Mergin')
    
    histList, titleList  = getListOfHistograms(InputDir,sampleDict)
    #logger.debug(f'List of histograms : {histList}')

    labelStyleDict = configDict.get("Groups")
    labelsToPlot   = list(labelStyleDict.keys())

    # print yield tables only
    if args.yieldOnly:
        yieldTables = getYields(sampleDict, InputDir, lumi)
        with open(os.path.join(plotDir,'YieldTables_'+str(args.yaml).split('.')[0]+'.txt'), 'w') as yieldf:
            for table in yieldTables:
                print(table)
                yieldf.write(str(table)+'\n')
    else:
        # To save the legends only
        canvOnly, legOnly = getCanvas("canv_legOnly", 
                                      [0.02,0.02,0.98,0.97], 
                                      [0.02,0.02,0.98,0.02], 
                                      [800,100], 
                                      legendTxtSize=0.3, 
                                      logy=False)
        
        # Loop over all the histos and samples to produce all the normalised plots
        outRoot  = ROOT.TFile(os.path.join(plotDir,'plots_'+str(args.yaml).split('.')[0]+'.root'),'RECREATE')
        for hIdx, hName in enumerate(histList):
            #logger.info(f'Histogram :...::.. {hName}')
            labelHistDict = defaultdict(list)
            for sampleKey, sampleValDict in sampleDict.items():
                #logger.debug(f'Sample Name : {sampleKey}')
                filepath = os.path.join(InputDir,sampleValDict.get('file'))
                if not os.path.isfile(filepath):
                    print(f'W A R N I N G :: -----> {infile} not found !!!')
                    continue
                label  = sampleValDict.get('label')
                rfile  = ROOT.TFile(filepath,'READ')
                if rfile.IsZombie():
                    logger.warning(f'{filepath} is a Zombie !!!')
                    continue
                hist = copy.deepcopy(rfile.Get(hName))
                hist.SetDirectory(0)
                rfile.Close()
                #logger.debug(f'{label} ---> {hist} --> Integral : {hist.Integral()}')
                if label in labelsToPlot:
                    labelHistDict[label].append(hist)
            #logger.debug(f'Dictionary >>---> key : Label, value : List of histograms :: \n{labelHistDict}')    

            canv1, legend = getCanvas("canv_"+hName, legendPos, mergin, canvSize, logy=False)
            canv2, _ = getCanvas("canv_log"+hName, legendPos, mergin, canvSize, logy=True)
            draw_opts = "hist"
            stack = ROOT.THStack()
            for pIdx, (lab, plotList) in enumerate(labelHistDict.items()):
                if len(plotList) == 0:
                    logger.critical(f'No histograms found under label : {lab}')
                    sys.exit()
                linecolor = labelStyleDict.get(lab).get('linecolor')
                linestyle = labelStyleDict.get(lab).get('linestyle')
                linewidth = labelStyleDict.get(lab).get('linewidth')
                histo = copy.deepcopy(getHistPerLabel(plotList))
                histo.SetName(hName+'_'+lab.replace('+',''))
                histo.SetLineWidth(linewidth)
                histo.SetLineColor(linecolor)
                histo.SetLineStyle(linestyle)
                outRoot.cd()
                histo.Write()
                #logger.debug(f'{pIdx} -> {lab} -> {histo} -> {histo.Integral()}')
                stack.Add(histo)
                legend.AddEntry(histo, str(lab), "l")
                if hIdx == 0:
                    canvOnly.cd()
                    legOnly.AddEntry(histo, str(lab), "l")
                

            canv1.cd()
            stack.Draw("nostack" + draw_opts) 
            #legend.Draw()
            canv1.Update()
            #stack.GetYaxis().SetTitle("#frac{1}{N} (#frac{dN}{dx})")
            stack.GetYaxis().SetTitle("1/N dN/dx")
            stack.GetYaxis().SetTitleSize(0.05)
            stack.GetYaxis().SetTitleOffset(1.11)
            stack.GetYaxis().CenterTitle(True)
            stack.GetYaxis().SetLabelSize(0.042)
            stack.GetYaxis().SetTickLength(0.02)
            stack.GetXaxis().SetTitle(titleList[hIdx])
            stack.GetXaxis().SetTitleSize(0.06)
            stack.GetXaxis().SetTitleOffset(0.87)
            stack.GetXaxis().CenterTitle(True)
            stack.GetXaxis().SetLabelSize(0.044)
            stack.GetXaxis().SetTickLength(0.02)
            canv1.Modified()
            canv1.SaveAs(os.path.join(plotDir,hName+'.'+plotExt))

            canv2.cd()
            stack.Draw("nostack" + draw_opts)
            #legend.Draw()
            canv2.Update()
            #stack.GetYaxis().SetTitle("#frac{1}{N} (#frac{dN}{dx}) [log scale]")
            stack.GetYaxis().SetTitle("1/N dN/dx [log scale]")
            stack.GetYaxis().SetTitleSize(0.05)
            stack.GetYaxis().SetTitleOffset(1.11)
            stack.GetYaxis().CenterTitle(True)
            stack.GetYaxis().SetLabelSize(0.042)
            stack.GetYaxis().SetTickLength(0.02)
            stack.GetXaxis().SetTitle(titleList[hIdx])
            stack.GetXaxis().SetTitleSize(0.06)
            stack.GetXaxis().SetTitleOffset(0.87)
            stack.GetXaxis().CenterTitle(True)
            stack.GetXaxis().SetLabelSize(0.044)
            stack.GetXaxis().SetTickLength(0.02)
            canv2.Modified()
            canv2.SaveAs(os.path.join(plotDir,hName+'_logy.'+plotExt))

        outRoot.Close()
        canvOnly.cd()
        legOnly.Draw()
        logger.info("Saving legend separately ... ")
        canvOnly.SaveAs(os.path.join(plotDir, "Legend."+plotExt))

if __name__ == "__main__":
    main()
