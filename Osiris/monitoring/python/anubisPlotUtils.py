import sys, os
import datetime
import json
import struct
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.image as image
import mplhep as hep
hep.style.use([hep.style.ATLAS])
from PIL import Image
from anubisDataUtils import *

def heatFromROOTFile(dataFile, time=1800, name="proANUBISMonitorHeatMap"):
    #Plots heat maps from triggered data, showing the hit rate in each rpc channel. 2D plots designed to replicate RPC layout and channel counting direction.
    dataDict = importFromROOTFile(dataFile)
    thisHitData = {}
    addresses = ['ee00','ee01','ee02','ee03','ee04']
    for tdc in range(5):
        thisHitData[addresses[tdc]] = countChannels(dataDict["data"][tdc])
    makeHitMaps(thisHitData,dataDict["endTime"], name,False,unit='hz',time=time)

def heatFromRawfile(dataFile, time=1800, name="proANUBISMonitorHeatMap"):
    #Plots heat maps from triggered data, showing the hit rate in each rpc channel. 2D plots designed to replicate RPC layout and channel counting direction.
    dataDict = importFromRawfile(dataFile)
    thisHitData = {}
    addresses = ['ee00','ee01','ee02','ee03','ee04']
    for tdc in range(5):
        thisHitData[addresses[tdc]] = countChannels(dataDict["data"][tdc])
    makeHitMaps(thisHitData,dataDict["endTimeStr"], name,False,unit='hz',time=time)    

def plotLumiEventRate(period, name="ATLAS_vs_proANUBIS_Lumi", outputDir=os.environ["MON_PLOT_DIR"]):

        # For a given time period:
        #   - Get all the datafiles in that range 
        #       * If loading files that fail, get the filename and note it as corrupted
        #   - Get the ATLAS Lumi in that range
        #       * Not sensible to keep a local copy of all ATLAS' lumi measurements.
        #       * Can we access the lumi info remotely as needed, without downloading particular periods?
        return 0
        

def plotTDCreads(datafile, name="TDCreads", outputDir=os.environ["MON_PLOT_DIR"]):
    dataDict = importDataFile(dataFile)

    if "wordCounts" not in dataDict.keys():
        print(f"{datafile} does not support plotTDCreads currently, requires 'wordCounts', try another file")
        return ""

    fig, ax = plt.subplots(1, 1, figsize=(12,8))

    ax.hist(dataDict["wordCounts"], bins=100, fill=False,histtype='step', linewidth=2.0)
    ax.set_yscale('log')
    ax.set_ylabel("TDC Reads")
    ax.set_xlabel("32-bit Words")
    plotName= f"{outputDir}/{name.strip(' ')}.png"
    plt.savefig(plotName)
    plt.close()
    return plotName 

def plotMaxRates(datafiles, time=1800, name="maxRates", splitTDC=False, outputDir=os.environ["MON_PLOT_DIR"]):
    dates= []
    maxRates= {"Total": []}
    suffix="" if not splitTDC else "_splitTDC" 
    for infile in args.inFile:
        
        tempMaxRates = maxRateFromFile(infile, time=time)
        
        dates.append(date)
        #TODO: If date appears multiple times only take the first one and take max rate of the different cases
        if splitTDC: 
            for tdc, rate in tempMaxRates.items():
                if tdc in maxRates.keys():
                    maxRates[tdc].append(rate)
                else:
                    maxRates[tdc]=rate
        else:
            maxRates["Total"].append(tempMaxRates["Total"])

    #TODO Break down the max rates by each TDC using splitTDC

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.step(dates,maxRates["Total"], where='post', linewidth=3)
    ax.set_yscale('log')
    ax.set_xlabel('Date')
    ax.set_ylabel('ANUBIS_DAQ Max Readout Rate (Hz)')
    #im = image.imread('/eos/user/m/mireveri/anubis/ANUBISLogo.png')
    im = image.imread(os.environ["ANUBIS_LOGO"])
    xScale = abs(max(dates)-min(dates))
    yScale = abs(max(maxRates["Total"])-min(maxRates["Total"]))
    ax.imshow(im, aspect='auto', extent=(0.1*xScale, 0.2*xScale, 0.7*yScale, 0.9*yScale), zorder=1)
    ax.set_xlim([min(dates),max(dates[-1])])
    ax.set_ylim(1,1.2*max(maxRates))
    plotName= f"{outputDir}/{name.strip(' ')}{suffix}.png"
    plt.savefig(plotName)
    plt.close()
    return plotName 

def plotPhi(array, name, zrange = [0.01,200], unit='khz', time=60., outputDir=os.environ["MON_PLOT_DIR"]):
    fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
    #lab = hep.atlas.label(com=False,data=True, label="Internal")
    #lab[2].set_text("")
    #im = image.imread('/eos/user/m/mireveri/anubis/ANUBISLogo.png')
    im = image.imread(os.environ["ANUBIS_LOGO"])
    phichannels = [x-0.5 for x in range(65)]
    phiHist = ((np.array([array])/time).transpose(),np.array(phichannels),np.array([0,1]))
    thisHist = hep.hist2dplot(phiHist,norm=colors.LogNorm(zrange[0],zrange[1]))
    thisHist.cbar.set_label('Event Rate ('+unit+')', rotation=270, y=0.4)
    plt.xlabel("Channel")
    plt.ylabel(" ")
    plt.title(name)
    fig.tight_layout()
    ax.get_yaxis().set_visible(False)
    ax.imshow(im, aspect='auto', extent=(1, 3, .88, .93), zorder=1)
    plotName= f"{outputDir}/{name.replace(' ','')}.png"
    #plt.savefig('/eos/user/m/mireveri/anubis/'+name.strip(" ")+".png")
    plt.savefig(plotName)
    plt.close()
    return plotName 

def plotEta(array, name, zrange = [0.01,200], unit='khz', time=60., outputDir=os.environ["MON_PLOT_DIR"]):
    fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
    #lab = hep.atlas.label(com=False,data=True, label="Internal")
    #lab[2].set_text("")
    #im = image.imread('/eos/user/m/mireveri/anubis/ANUBISLogo.png')
    im = image.imread(os.environ["ANUBIS_LOGO"])
    etachannels = [x-0.5 for x in range(33)]
    etaHist = (np.array([array])/time,np.array([0,1]),np.array(etachannels))
    thisHist = hep.hist2dplot(etaHist,norm=colors.LogNorm(zrange[0],zrange[1]))
    thisHist.cbar.set_label('Event Rate ('+unit+')', rotation=270, y=0.4)
    plt.ylim(31.5,-0.5)
    plt.ylabel("Channel")
    plt.xlabel(" ")
    plt.title(name)
    fig.tight_layout()
    ax.get_xaxis().set_visible(False)
    ax.imshow(im, aspect='auto', extent=(0.03, 0.06, 3.2, 1.6), zorder=1)
    plotName=f"{outputDir}/{name.replace(' ','')}.png"
    #plt.savefig('/eos/user/m/mireveri/anubis/'+name.strip(" ")+".png")
    plt.savefig(plotName)
    plt.close()
    return plotName 

def combinePlots(plots,imname, outputDir=os.environ["MON_PLOT_DIR"]):
    images = [Image.open(x) for x in plots]
    widths, heights = zip(*(i.size for i in images))

    total_width = int(2*widths[0])
    if(len(plots)|2>0):
        max_height = int((sum(heights)+heights[0])/2)
    else:
        max_height = int(sum(heights)/2)

    new_im = Image.new('RGB', (total_width, max_height))

    x_offset = 0
    y_offset = 0
    even = True
    for im in images:
        if even:
            new_im.paste(im, (x_offset,y_offset))
            x_offset += im.size[0]
            even = False
        else:
            new_im.paste(im,(x_offset,y_offset))
            x_offset = 0
            y_offset += im.size[1]
            even = True

    new_im.save(f"{outputDir}/{imname.strip(' ')}.pdf")
    
def makeHitMaps(filenames, plotTitle, imname, useJson=True, unit='khz', time=60., outputDir=os.environ["MON_PLOT_DIR"]):
    if(useJson):
        hitData = {}
    
        addresses = ['ee00','ee01','ee02','ee03','ee04']
        for fname in filenames:
            thisFile = open(fname)
            jsonData = json.load(thisFile)
            for addr in addresses:
                try:
                    hitData[addr]=jsonData['Summary']['TDCs'][addr]['nHits']
                except KeyError:
                    continue
        for addr in addresses:
            if addr not in hitData.keys():
                hitData[addr] = [0 for x in range(128)]
            else:
                for idx, hit in enumerate(hitData[addr]):
                    hitData[addr][idx]= hit/1000. #Divide by 1000 to convert to khz
    else:
        hitData = filenames
    tripEtaLow = hitData['ee00'][0:32]
    tripPhiLow = hitData['ee00'][32:96]
    tripEtaMid = hitData['ee00'][96:128]
    tripPhiMid = hitData['ee01'][0:64]
    tripEtaTop = hitData['ee01'][64:96]
    tripPhiTop = hitData['ee01'][96:128]+hitData['ee02'][0:32]
    singEta = hitData['ee02'][32:64]
    singPhi = hitData['ee02'][64:128]
    doubEtaLow = hitData['ee03'][0:32]
    doubPhiLow = hitData['ee03'][32:96]
    doubEtaTop = hitData['ee03'][96:128]
    doubPhiTop = hitData['ee04'][0:64]
    imageArr = []
    imageArr.append(plotPhi(tripPhiLow,"Phi Triplet Low "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    imageArr.append(plotEta(tripEtaLow,"Eta Triplet Low "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    imageArr.append(plotPhi(tripPhiMid,"Phi Triplet Mid "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    imageArr.append(plotEta(tripEtaMid,"Eta Triplet Mid "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    imageArr.append(plotPhi(tripPhiTop,"Phi Triplet Top "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    imageArr.append(plotEta(tripEtaTop,"Eta Triplet Top "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    imageArr.append(plotPhi(singPhi,"Phi Singlet "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    imageArr.append(plotEta(singEta,"Eta Singlet "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    imageArr.append(plotPhi(doubPhiLow,"Phi Doublet Low "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    imageArr.append(plotEta(doubEtaLow,"Eta Doublet Low "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    imageArr.append(plotPhi(doubPhiTop,"Phi Doublet Top "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    imageArr.append(plotEta(doubEtaTop,"Eta Doublet Top "+plotTitle, unit=unit, time=time, outputDir=outputDir))
    combinePlots(imageArr,imname)
    for im in imageArr:
        os.remove(im)

def makeSingleLayer(data, name, outputDir=os.environ["MON_PLOT_DIR"]):
    #Heatmap plot of one RPC layer. Takes already-split heat map, used by event display
    fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
    channels= [x-0.5 for x in range(len(data)+1)]
    if(len(data)==32):
        histArr = (np.array([data]),np.array([0,1]),np.array(channels))
    else:
        histArr = ((np.array([data])).transpose(),np.array(channels),np.array([0,1]))
    thisHist = hep.hist2dplot(histArr,norm=colors.LogNorm(0.1,2))
    thisHist.cbar.remove()
    if(len(data)==32):
        plt.ylim(len(data)-0.5,-0.5)
    plt.ylabel(" ")
    plt.xlabel(" ")
    #plt.title(name)
    
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plotName = f"{outputDir}/{name}.png"
    #plt.savefig("/eos/user/m/mireveri/anubis/"+name+".png")
    plt.savefig(plotName)
    plt.close()
    return plotName 

def stackAndConvert(images, name="testDisplay", outputDir=os.environ["MON_PLOT_DIR"]):
    #PIL hacking to distort images and put them together to make a primitive replication of the detector
    img = Image.open(images[0])
    total_width = 3*img.size[0]
    max_height = int(4*img.size[1])
    new_im = Image.new('RGB', (total_width, max_height))
    newData = new_im.load()
    x_offset = 0
    y_offset = 6*int(max_height/8.)
    for y in range(max_height):
        for x in range(total_width):
            #newData forms the background of the image, need to set it to all-white to start. Probably some better way to do this?
            newData[x, y] = (255, 255, 255)
    for idx, image in enumerate(images):
        img = Image.open(image)
        img = img.convert("RGBA")
        temp_im = Image.new('RGBA', (3*img.size[0], img.size[1]))
        temp_im.paste(img, (int(img.size[0]/2.),0))
        temp_im = temp_im.transform(temp_im.size, Image.AFFINE, (0.5, 1., 0, 0, 1, 0))
        pixdata = temp_im.load()
        width, height = temp_im.size
        for y in range(height):
            for x in range(width):
                if pixdata[x, y] == (255, 255, 255, 255):
                    #Manually make any white pixel transparent so that they can stack together nicely.
                    pixdata[x, y] = (255, 255, 255, 0)
        new_im.paste(temp_im, (0, y_offset), temp_im)
        y_offset = y_offset-int(max_height/28.)
        if idx == 5 or idx==7:
            #Counts from the bottom up, want bigger gaps between the different chambers
            y_offset = y_offset-5*int(max_height/28.)                   
    #new_im.save("/eos/user/m/mireveri/anubis/"+name.strip(" ")+".png", "PNG")
    new_im.save(f"{outputDir}/{name.strip(' ')}.png", "PNG")
    
def makeEventDisplay(eventData,name, outputDir=os.environ["MON_PLOT_DIR"]):
    #Expects a single event, divided as [tdc0,tdc2,...,tdc4]
    countOne = countChannels([eventData[0]])
    countTwo = countChannels([eventData[1]])
    countThree = countChannels([eventData[2]])
    countFour = countChannels([eventData[3]])
    countFive = countChannels([eventData[4]])
    singEventPlots = []
    singEventPlots.append(makeSingleLayer(countOne[0:32],"Eta Triplet Low, Three Coincidences Required", outputDir=outputDir))
    singEventPlots.append(makeSingleLayer(countOne[32:96],"Phi Triplet Low, Three Coincidences Required", outputDir=outputDir))
    singEventPlots.append(makeSingleLayer(countOne[96:128],"Eta Triplet Mid, Three Coincidences Required", outputDir=outputDir))
    singEventPlots.append(makeSingleLayer(countTwo[0:64],"Phi Triplet Mid, Three Coincidences Required", outputDir=outputDir))
    singEventPlots.append(makeSingleLayer(countTwo[64:96],"Eta Triplet Top, Three Coincidences Required", outputDir=outputDir))
    singEventPlots.append(makeSingleLayer(countTwo[96:128]+countThree[0:32],"Phi Triplet Top, Three Coincidences Required", outputDir=outputDir))
    singEventPlots.append(makeSingleLayer(countThree[32:64],"Eta Singlet, Three Coincidences Required", outputDir=outputDir))
    singEventPlots.append(makeSingleLayer(countThree[64:128],"Phi Singlet, Three Coincidences Required", outputDir=outputDir))
    singEventPlots.append(makeSingleLayer(countFour[0:32],"Eta Doublet Low, Three Coincidences Required", outputDir=outputDir))
    singEventPlots.append(makeSingleLayer(countFour[32:96],"Phi Doublet Low, Three Coincidences Required", outputDir=outputDir))
    singEventPlots.append(makeSingleLayer(countFour[96:128],"Eta Doublet Top, Three Coincidences Required", outputDir=outputDir))
    singEventPlots.append(makeSingleLayer(countFive[0:64],"Phi Doublet Top, Three Coincidences Required", outputDir=outputDir))
    stackAndConvert(singEventPlots,name)
    for plot in singEventPlots:
        #Remove all the temporary plots. There's probably a better way to move pil images around than making and deleting .png files.
        os.remove(plot)
    return

def GetEvent(eventData, num):
    return [eventData[0][num],eventData[1][num],eventData[2][num],eventData[3][num],eventData[4][num]] 
