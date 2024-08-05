import sys, os
import datetime
import json
import struct
import ROOT
import numpy as np
from anubisPlotUtils import *
from glob import glob
import re

#===========================#
#   Data Importing Scripts  #
#                           #
#===========================#

def importDatafile(filename):
    if "txt" in filename.split(".")[-1]:
        return importFromTextFile(filename)
    elif "h5" in filename.split(".")[-1]:
        return importFromHDF5File(filename)
    elif "root" in filename.split(".")[-1]:
        return importFromROOTFile(filename)
    elif "raw" in filename.split(".")[-1]:
        return importFromRawFile(filename)
    else:
        print(f"File type ({filename.split('.')[-1]}) not recognized. Expect .txt, .root, .raw, or .h5 input file.")

def importFromTextFile(filename):
    inputText = open(filename)
    thisEvent = []
    data = [[] for tdc in range(5)]
    tdc=0
    for line in inputText:
        if "Header" in line:
            thisEvent = []
            tdc = int(line.split(" ")[1].strip(","))
        elif "Data" in line:
            thisEvent.append(int("0x"+line.split(" ")[2].strip("."),0))
        elif "EOB" in line:
            data[tdc].append(thisEvent)
    return {"data": data}

def importFromHDF5File(filename):
    inputHDF5 = h5py.File(filename)
    thisEvent = []
    data = [[] for tdc in range(5)]
    tdc=0
    for event in inputHDF5['data']:
        tdc = event[0]-60928
        thisEvent = []
        for hit in event[2]:
            thisEvent.append(hit)
        data[tdc].append(thisEvent)
    
    return {"data": data}

def importFromROOTFile(filename, style = "New"):
    inputRoot = ROOT.TFile(filename, "READ")
    eventTree = inputRoot.Get("Events")    
    data = [[] for tdc in range(5)]
    nevts = [0 for tdc in range(5)]
    eventTree.GetEvent(0)
    if "New" not in style:
        lastTDC = eventTree.TDC
        for event in range(eventTree.GetEntries()):
            eventTree.GetEvent(event)
            tdc = eventTree.TDC
            if tdc is not lastTDC:
                if nevts[lastTDC]!= max(nevts):
                    print("Missed trigger on TDC", lastTDC, "total event number",event)
                    data[lastTDC].append([[],-1])
                    nevts[lastTDC] = nevts[lastTDC]+1
                lastTDC=tdc
            thisEvent = []
            hits = eventTree.Words
            for hit in hits:
                thisEvent.append(hit)
            data[tdc].append([thisEvent,eventTree.TimeStamp.AsDouble()])
            nevts[tdc] = nevts[tdc]+1
    else:
        headerMask = 0x400000
        EOBMask = 0x200000
        for event in range(eventTree.GetEntries()/2.,eventTree.GetEntries()/2.+10):
            eventTree.GetEntry(event)
            thisEvent = []
            thisEventNum = 0
            qualFlag = 0
            ready=False
            for word in eventTree.Words:
                if word&headerMask:
                    if not ready:
                        ready=True
                        thisEventNum = word&0xfff
                        qualFlag=0
                        thisEvent = []
                    else:
                        if((word&0xfff)==thisEventNum+1):
                            qualFlag = qualFlag|0x2
                            data[eventTree.TDC].append([thisEvent,eventTree.TimeStamp.AsDouble(), qualFlag])
                            thisEvent = []
                        else:
                            qualFlag = qualFlag|0x1
                elif word&EOBMask:
                    if ready:
                        data[eventTree.TDC].append([thisEvent,eventTree.TimeStamp.AsDouble(), qualFlag])
                        qualFlag = 0
                        thisEvent = []
                        ready=False
                else:
                    if not ready:
                        ready=True
                        tdc=eventTree.TDC
                        thisTime=eventTree.TimeStamp.AsDouble()
                        thisEvent.append(word)
                        qualFlag = qualFlag|0x4
                    else:
                        thisEvent.append(word)
            if ready:
                qualFlag = qualFlag|0x2
                data[eventTree.TDC].append([thisEvent,eventTree.TimeStamp.AsDouble(), qualFlag])
    endTime = inputRoot.Get("endTime").AsString()
    
    return {"data": data, "endTime": endTime}

def importFromRawfile(filename):
    outdata = [[] for tdc in range(6)]
    wordCounts = []
    readTimes = []
    with open(filename,'rb') as data:
        headerMask = 0x400000
        EOBMask = 0x200000
        dat = data.read()
        words = struct.iter_unpack("I",dat)
        state=0
        ready = False
        nwords = 0
        thisTDC = 0
        startTimeSec = 0
        endTimeSec = 0
        startTimeNS = 0
        endTimeNS = 0
        thisTimeSec = 0
        thisTimeNS = 0
        qualFlag = 0
        lastWasOne = False
        thisEvent = []
        thisEventNum=0
        nEvents = 0
        for word in words:
            thisWord = word[0]
            if state==0:
                startTimeSec = thisWord
                state=1
            elif state==1:
                startTimeNS = thisWord
                state=2
            elif state==2:
                thisTimeSec = thisWord
                state=3
            elif state==3:
                thisTimeNS = thisWord
                state=4
            elif state==4:
                thisTDC = (thisWord>>24)
                nwords = thisWord&0xffffff
                if(thisTDC<5):
                    timestamp = thisTimeSec+thisTimeNS/1000000000.
                    readTimes.append(datetime.datetime.fromtimestamp(timestamp))
                    wordCounts.append(nwords)
                    if(nwords==1):
                        lastWasOne=True
                        #print(hex(thisWord))
                state=5
            elif state==5:
                if(lastWasOne):
                    #print(hex(thisWord))
                    lastWasOne=False
                if thisWord&headerMask:
                    if not ready:
                        ready=True
                        thisEventNum = thisWord&0xfff
                        qualFlag=0
                        thisEvent = []
                    else:
                        if((thisWord&0xfff)==thisEventNum+1):
                            qualFlag = qualFlag|0x2
                            #outdata[thisTDC].append(thisEvent)
                            if(thisTDC==0):
                                nEvents = nEvents+1
                            thisEvent = []
                        else:
                            qualFlag = qualFlag|0x1
                elif thisWord&EOBMask:
                    if ready:
                        outdata[thisTDC].append(thisEvent)
                        if(thisTDC==0):
                            nEvents = nEvents+1
                        qualFlag = 0
                        thisEvent = []
                        ready=False
                else:
                    if not ready:
                        ready=True
                        thisEvent.append(thisWord)
                        qualFlag = qualFlag|0x4
                    else:
                        thisEvent.append(thisWord)          
                    if(thisTDC==5): #Runs in continuous mode without Header or EOB
                        outdata[thisTDC].append([thisWord])

                nwords = nwords-1
                if nwords==0:
                    state=2
    print(nEvents)
    print(datetime.datetime.fromtimestamp(startTimeSec).strftime("%m_%d_%Y_%H:%M:%S"))
    print(datetime.datetime.fromtimestamp(thisTimeSec).strftime("%m_%d_%Y_%H:%M:%S"))
    outputDict={"data": outdata,
                "startTimeStr": datetime.datetime.fromtimestamp(startTimeSec).strftime("%m_%d_%Y_%H:%M:%S"),
                "endTimeStr": datetime.datetime.fromtimestamp(thisTimeSec).strftime("%m_%d_%Y_%H:%M:%S"),
                "startTime": datetime.datetime.fromtimestamp(startTimeSec),
                "endTime": datetime.datetime.fromtimestamp(thisTimeSec),
                "wordCount": wordCounts,
                "readTimes": readTimes,
                }

    return outputDict

#---------------------------------------------------------------------------------------------------------------

#===========================#
#   Data Handling Scripts   #
#                           #
#===========================#

def getStartAndEndTimesFromRawfile(filename):
    with open(filename,'rb') as data:
        startTimeSec = 0
        endTimeSec = 0
        startTimeNS = 0
        endTimeNS = 0
        dat = data.read(8)
        words = struct.iter_unpack("I",dat)
        for idx, word in enumerate(words):
            thisWord = word[0]
            startTimeSec = startTimeNS
            startTimeNS = thisWord
            if idx>0:
                break
        data.seek(-8, os.SEEK_END)
        dat = data.read()
        words = struct.iter_unpack("I",dat)
        for idx, word in enumerate(words):
            thisWord = word[0]
            endTimeSec = endTimeNS
            endTimeNS = thisWord
        startTime = datetime.datetime.fromtimestamp(startTimeSec+startTimeNS/1000000000.)
        endTime = datetime.datetime.fromtimestamp(endTimeSec+endTimeNS/1000000000.)

    outputDict={"startTime": startTime,
                "endTime": endTime
                }

    return outputDict

def getDuration(filename):
    if ".raw" in filename:
        timeDict = getStartAndEndTimesFromRawfile(filename)
        duration = (timeDict['endTime'] - timeDict['startTime']).total_seconds()
    else:
        raise Exception(f"Not implmented for this file type: .{filename.split('.')[-1]}")

    return duration

def printTotalPerChannel(filename, doRate = False):
    if doRate:
        print(f"Rates:")
        unit="Hz"
    else:
        print(f"Counts:")
        unit=""

    dataDict = importFromRawfile(filename)
    thisHitData = {}
    for tdc in range(6):
        print(f"TDC {tdc}")
        thisHitData = countChannels(dataDict["data"][tdc])

        for channel in range(len(thisHitData)):
            countOrRate = thisHitData[channel]
            if doRate: 
                countOrRate/=duration

            outString=f"{countOrRate}" if float(countOrRate)==0 else '\033[91m'+f"{countOrRate}"+'\033[0m'
            print(f"Ch{channel}: {outString}{unit}, ", end="")
            if (int(channel)+1)%8==0:
                print("")

def countChannels(events):
    #Expects events from one TDC, counts how many hits each channel has within the event list
    chanCounts = [0 for x in range(128)]
    for event in events:
        for word in event:
            try:
                chanCounts[(word>>24)&0x7f]+=1
            except:
                print(word>>24)
    return chanCounts

def maxRateFromFile(dataFile, time=1800):
    dataDict = importDataFile(dataFile)
    maxRates = {'Total': -1}
    addresses = ['ee00','ee01','ee02','ee03','ee04']
    
    for tdc in range(5):
        rates = countChannels(dataDict["data"][tdc])/time
        maxRates[addresses[tdc]] = max(rates)
        if max(rates) > maxRates['Total']:
            maxRates['Total'] = max(rates)

    return maxRates

# Get list of timestamps in file:
#   Assumes file format: *_YYMMDD_hhmm.*
def getTimestamp(file):
    timeStamp = {"timeStamp": "", "date": "", "time": ""}

    timeStamp["timeStamp"]= next(re.finditer(r"\d+_\d+", file.split("/")[-1], re.MULTILINE)).group(0)
    timeStamp["date"]= timeStamp["timeStamp"].split("_")[0]
    timeStamp["time"]= timeStamp["timeStamp"].split("_")[1]

    return timeStamp

# Get the most recent Data file
def getMostRecentData(baseDir="/eos/atlas/atlascerngroupdisk/proj-anubis/proANUBIS/data/"):
    currentMonth=datetime.datetime.today().strftime('%y_%m')
    currentDate=datetime.datetime.today().strftime('%y%m%d')

    recentDataList=[]
    for fileFormats in ["raw","ROOT"]:
        recentDataList.extend(glob(baseDir+f"{currentMonth}/{currentDate}/*.{fileFormats}"))

    if len(recentDataList)==0:
        #Look for other days in that month
        print(f"Looking for files in {currentMonth}")
        for fileFormats in ["raw","ROOT"]:
            recentDataList.extend(glob(baseDir+f"{currentMonth}/*/*.{fileFormats}"))

    if len(recentDataList)==0:
        #Look in the previous month too
        previousMonth = (datetime.datetime.today().replace(day=1) - datetime.timedelta(days=1))
        print(f"Looking for files in {previousMonth.strftime('%y_%m')}")
        for fileFormats in ["raw","ROOT"]:
            recentDataList.extend(glob(baseDir+f"{previousMonth.strftime('%y_%m')}/*/*.{fileFormats}"))

    mostRecent="-1"
    for recentData in recentDataList:
        if mostRecent == "-1":
            mostRecent = recentData
        else:
            tempTimeStamp= getTimestamp(mostRecent)
            timeStamp= getTimestamp(recentData)

            if timeStamp["time"] > tempTimeStamp["time"]:
                mostRecent = recentData

    print(f"Found most recent file: {mostRecent}")

    return mostRecent
