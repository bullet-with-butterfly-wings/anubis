import argparse
import struct
import re

def countCorruption(filename, ntdcs, doPrint):
    missedEvents = [[] for tdc in range(ntdcs)]
    presentEvents = [[] for tdc in range(ntdcs)]
    lastWord = [-1 for tdc in range(ntdcs)]
    extraHeads = [0 for tdc in range(ntdcs)]
    extraTails = [0 for tdc in range(ntdcs)]
    nHits = [0 for tdc in range(ntdcs)]
    corruptedHeaders = [0 for tdc in range(ntdcs)]
    with open(filename,'rb') as data:
        headerMask = 0x400000
        EOBMask = 0x200000
        dat = data.read()
        words = struct.iter_unpack("I",dat)
        state=0
        thisTDC = -1
        startTimeSec = 0
        endTimeSec = 0
        startTimeNS = 0
        endTimeNS = 0
        thisTimeSec = 0
        thisTimeNS = 0
        lastWasHead = False
        lastWasTails = False
        for idx, word in enumerate(words):
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
                if doPrint:
                    print("Block Read Header from TDC",thisTDC,", nwords =",nwords);
                firstEvt = -1
                lastEvt = 0
                nHeads = 0
                lastWasHead = False
                lastWasTails = False
                state=5
                if(thisTDC==5):
                    state=6
                if(thisTDC)>6:
                    print(idx)
                    break
            elif state==5:
                if doPrint:
                    if thisWord&headerMask:
                        print("TDC",thisTDC,", Header:",hex(thisWord))
                    elif thisWord&EOBMask:
                        print("TDC",thisTDC,", EOB:",hex(thisWord))
                    else:
                        print("TDC",thisTDC,", Data word:",hex(thisWord))

                if thisWord&headerMask:
                    thisEventNum = thisWord&0x3ff
                    if firstEvt<0:
                        firstEvt = thisEventNum
                    lastEvt = thisEventNum
                    if lastWasHead:
                        extraHeads[thisTDC] = extraHeads[thisTDC]+1
                    else:
                        nHeads = nHeads+1
                        lastWasHead=True
                    lastWasTails = False
                else:
                    lastWasHead = False
                    if thisWord&EOBMask:
                        if lastWasTails:
                            extraTails[thisTDC] = extraTails[thisTDC]+1
                        else: 
                            lastWasTails=True
                    else:
                        nHits[thisTDC]+=1
                        lastWasTails = False
                nwords = nwords-1
                if nwords==0:
                    if lastWord[thisTDC]>-0.5 and firstEvt>-0.5:
                        evtDiff = firstEvt - lastWord[thisTDC]
                        if evtDiff < 1:
                            evtDiff = evtDiff+0x3ff+1

                        if evtDiff != 0:
                            missedEvents[thisTDC].append(evtDiff - 1)
                        else:
                            corruptedHeaders[thisTDC]+=1

                    lastWord[thisTDC] = lastEvt
                    presentEvents[thisTDC].append(nHeads)

                    state=2
            elif state==6:
                nwords = nwords-1
                if nwords==0:
                    state=2
    outDict = {'missedEvents':missedEvents,
			   'presentEvents':presentEvents,
			   'extraHeads':extraHeads,
			   'extraTails':extraTails,
			   'nHits': nHits,
               'corruptedHeaders': corruptedHeaders,
               }
    return outDict

def checkAdditionals(file, addHits):
    lastWord = -1
    firstEvt = -1
    lastEvt = -1
    extraHeads = 0  
    extraTails = 0

    nHeads = 0
    nTails = 0
    nHits=0
    nStatus={}

    eventOverflow = 0

    with open(file, "r") as f:
        lines = f.readlines()
        lastWasHead = False
        lastWasTails = False

        for line in lines:

            if "HEADER" in line:
                thisEventNum = int(line.split("#=")[-1])
                if firstEvt<0:
                    firstEvt = thisEventNum
                
                if lastEvt==1024:
                    eventOverflow+=1

                lastEvt = thisEventNum
                if lastWasHead:
                    extraHeads+=1
                else:
                    nHeads+=1
                    lastWasHead=True
                lastWasTails = False
            else:
                lastWasHead = False
                if "EOB" in line:
                    status = int(line.split("Status=")[-1])
                    if lastWasTails:
                        extraTails+=1
                    else:
                        lastWasTails=True
                else:
                    if "EDGE" not in line: 
                        nHits+=1
                    lastWasTails = False
        
        expectedEvents = (firstEvt + nHeads % 1024) - 1 
        if expectedEvents > 1024:
            expectedEvents = (expectedEvents % 1024) -1

        missedEvents=abs(expectedEvents - lastEvt)

    print(f"  There were {nHeads} Headers")
    print(f"  There were {nHits} Hits")
    print(f"  There were {extraHeads} extra Headers")
    print(f"  There were {extraTails} extra EOBs")
    print(f"  There were {missedEvents} missed events\n")
    if addHits:
        print(f"  Hits, Headers, extraHeaders, extraTails, missedEvents")
        print(f"  {nHits} {nHeads} {extraHeads} {extraTails} {missedEvents}\n")
    else:
        print(f"  Headers, extraHeaders, extraTails, missedEvents")
        print(f"  {nHeads} {extraHeads} {extraTails} {missedEvents}\n")


parser = argparse.ArgumentParser(description='Provide an input file')
parser.add_argument('--inDir', default="", help='Input directory')
parser.add_argument('-f', nargs="+", default=[], help='A list of files to read')
parser.add_argument('-n', default = 1, help='The number of tdcs to run')
parser.add_argument('--ID', type=int, default = -1, help='The ID of job to run')
parser.add_argument('--minID', type=int, default = -1, help='The minimum ID of job to run')
parser.add_argument('-p', action='store_true')
parser.add_argument('--addHits', action='store_true')

args = parser.parse_args()

from glob import glob

if args.inDir !="":
    files = glob(args.inDir+"/*")
else:
    files = args.f

for file in files:
    if "SRb" in file: 
        continue
    if ".txt" in file:
        if "ID" in file:
            ID = int(file.split("_ID")[-1].replace(".txt",""))
            if ID < args.minID:
                continue

            if ID != args.ID and args.ID!=-1:
                continue

        print(file)
        print("Run Summary:")
        checkAdditionals(file, args.addHits)
    else:
        ntdcs = args.n

        if "ID" in file:
            ID = int(file.split("_ID")[-1].replace(".raw",""))
            if ID < args.minID:
                continue

            if ID != args.ID and args.ID!=-1:
                continue

            if ( ID in [5,6,7,16,21,22,23,24,25,26,27,28,32,33,91,92,93,94,95] ) or (ID > 134 and ID < 160) :
                ntdcs = 2
            elif ID in [34,35,96,97,98,99,100,101]:
                ntdcs = 3
            elif ID in [36,37]:
                ntdcs = 4

        print(file)
        evtInfo = countCorruption(file, ntdcs, args.p)
        print("Run Summary:")
        for tdc in range(ntdcs):
            print("TDC",str(tdc))
            print("  Number of Hits:",evtInfo['nHits'][tdc],"")
            print("  Normal headers:",sum(evtInfo['presentEvents'][tdc]))
            print("  Extra headers:",evtInfo['extraHeads'][tdc])
            print("  Extra EOB:",evtInfo['extraTails'][tdc])
            print("  Missing Events:",sum(evtInfo['missedEvents'][tdc]))
            print("  Potentially Corrupted Headers:",evtInfo['corruptedHeaders'][tdc],"\n")

            if args.addHits:
                print(f"  Hits, Headers, extraHeaders, extraTails, missedEvents, corruptedHeaders")
                print(f"  {evtInfo['nHits'][tdc]} {sum(evtInfo['presentEvents'][tdc])} {evtInfo['extraHeads'][tdc]} {evtInfo['extraTails'][tdc]} {sum(evtInfo['missedEvents'][tdc])} {evtInfo['corruptedHeaders'][tdc]}\n")
            else:
                print(f"  Headers, extraHeaders, extraTails, missedEvents, corruptedHeaders")
                print(f"  {sum(evtInfo['presentEvents'][tdc])} {evtInfo['extraHeads'][tdc]} {evtInfo['extraTails'][tdc]} {sum(evtInfo['missedEvents'][tdc])} {evtInfo['corruptedHeaders'][tdc]}\n")

    print("=====================------------------------------====================\n")
    



