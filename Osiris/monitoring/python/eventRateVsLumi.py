from PIL import Image
import h5py
import anubisPlotUtils as anPlot
import json
import numpy as np
import os
import hist as hi
import ROOT
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.image as image
import mplhep as hep
import os
hep.style.use([hep.style.ATLAS])
import sys
import datetime
import csv
import matplotlib.dates as mdates
import struct
from pathlib import Path
from glob import glob

#baseDir="/eos/atlas/atlascerngroupdisk/proj-anubis/proANUBIS/data/24_04/"
baseDir="/eos/atlas/atlascerngroupdisk/proj-anubis/proANUBIS/data/24_05/"
#lumiFile=os.environ["MON_DIR"]+"/python/Luminosity_10e30-data-2024-04-19_10_56_18.csv"
lumiFile=os.environ["MON_DIR"]+"/python/Luminosity_10e30-data-2024-05-10_08_02_10.csv"

def getStartAndEndDate(dates):
    datetimeDates = [ datetime.datetime.strptime(date, "%y_%m_%d") for date in dates ] 
    startDate = (min(datetimeDates)).strftime("%y_%m_%d")
    endDate = (max(datetimeDates)).strftime("%y_%m_%d")

    return startDate, endDate

dates=[]
#for day in range(12,20):
#for day in range(3,13):
for day in range(3,11):
    #dates.append(f"24_04_{day}")
    if day < 10: 
        dates.append(f"24_05_0{day}")
    else:
        dates.append(f"24_05_{day}")

dataFiles=[]
for date in dates:
    tempFiles=glob(f"{baseDir}{date}/*")
    #if date == "24_04_12":
    if date == "24_05_03":
        for tF in tempFiles:
            if int(tF.split("_")[-1].split(".")[0]) > 1330:
                dataFiles.append(tF)
    else:
        dataFiles.extend(tempFiles)

xvals, yvals = [], []
failedFiles = []
for dF in dataFiles:
    print(dF)
    if ".root" in dF.lower():
        try:
            thisFile = ROOT.TFile(dF, "READ")
            start = datetime.datetime.fromtimestamp(thisFile.Get("startTime").AsDouble()-60*60)
            end = datetime.datetime.fromtimestamp(thisFile.Get("endTime").AsDouble()-60*60)
            if (start.year != 2024 or end.year != 2024):
                print("Invalid year in start or end time")
                print(f"Start Time: {start}")
                print(f"End Time: {end}")
                failedFiles.append(dF)
                continue
            xvals.append(start)
            yvals.append(1)
            xvals.append(end)
            yvals.append(0)
        except Exception as error:
            failedFiles.append(dF)
            print(f"An error occured: {type(error).__name__} -> {error}")
            continue
    elif ".raw" in dF.lower():
        try:
            #outDict= anPlot.importFromRawfile(dF)
            outDict= anPlot.getStartAndEndTimesFromRawfile(dF)
            if (outDict["startTime"].year != 2024 or outDict["endTime"].year != 2024):
                print("Invalid year in start or end time")
                print(f"Start Time: {outDict['startTime']}")
                print(f"End Time: {outDict['endTime']}")
                failedFiles.append(dF)
                continue
            xvals.append(outDict["startTime"])
            yvals.append(1)
            xvals.append(outDict["endTime"])
            yvals.append(0)
        except Exception as error:
            failedFiles.append(dF)
            print(f"An error occured: {type(error).__name__} -> {error}")
            continue

newLumis = []
counter = 10
with open(lumiFile, newline='') as csvfile:
    lumireader = csv.reader(csvfile, delimiter=',')
    for row in lumireader:
        if "Time" not in row[0]:
            newLumis.append((datetime.datetime.fromisoformat(row[0]),float(row[1])))
            if counter <10:
                print((datetime.datetime.fromisoformat(row[0]),float(row[1])))
                counter+=1

newLumis.sort(key = lambda x: x[0])


#-----------------#
#     Plotting    #
#-----------------#
startDate, endDate = getStartAndEndDate(dates)
outDir=""
outfileName=f"{outDir}proANUBISuptime-S{startDate}-E{endDate}.pdf"


fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
lumix = ax.twinx()
p1, = ax.step(xvals, yvals, 'k-', label='proANUBIS Uptime', where='post')
ax.set_xlabel("Date")
ax.set_ylabel("proANUBIS Data-taking State")
ax.fill_between(xvals,yvals,0,step="post",color="b")
#plt.xlim([splitTimes[1],splitTimes[48*3-12]])
p2, = lumix.plot([lumis[0] for lumis in newLumis], [lumis[1] for lumis in newLumis], 'r-', label='ATLAS Luminosity')
lumix.set_ylabel("ATLAS Luminosity / 10e30")
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
plt.gca().xaxis.set_major_locator(mdates.DayLocator())
#ax.legend(handles=[p1, p2])

#Adding the ANUBIS logo to the plot with scaled x length
im = image.imread(os.environ["ANUBIS_LOGO"])
xScale = max(xvals)-min(xvals)
print(min(xvals), max(xvals))
print(xScale, (0.82*xScale)+min(xvals), (0.9*xScale)+max(xvals))
ax.imshow(im, aspect='auto', extent=((0.82*xScale)+min(xvals), (0.9*xScale)+min(xvals), 1.05, 1.15), zorder=1)

ax.set_xlim(min(xvals),max(xvals))
ax.set_ylim([0,1.2])
ax.set_yticks([0, 1], ['Off', 'On'],rotation=0)  # Set text labels and properties.
#plt.show()
plt.savefig(outfileName)

print(f"Saved plot to: {outfileName}")

print("Failed Files:")
print(failedFiles)
