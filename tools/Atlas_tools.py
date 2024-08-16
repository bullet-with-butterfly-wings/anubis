import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
hep.style.use([hep.style.ATLAS])
import sys
import os
import matplotlib.colors as colors
import matplotlib.backends.backend_pdf
import hist as hi
import datetime
import importlib
import uproot

#sys.path.insert(1, 'C://Users//Peter//OneDrive - University of Cambridge//Desktop//summer2//Osiris Temp//processing//python')
sys.path.append(os.path.join(sys.path[0], 'Osiris Temp', 'processing', 'python'))
import rawFileReader
import Reconstruction_tools as RTools
import Timing_tools as TTools
importlib.reload(RTools)
importlib.reload(TTools)
from scipy.optimize import curve_fit

#scp -o "ProxyJump jd2052@gw.hep.phy.cam.ac.uk" jd2052@pcls.hep.phy.cam.ac.uk:/r04/atlas/revering/data/24_07/24_07_30/proAnubis_240730_2017.raw C:\Users\jony\Downloads

class BCR():
    def __init__(self, hitTime, event_time, bcr_count):
        self.hitTime = hitTime
        self.event_time = event_time
        self.absTime = event_time+datetime.timedelta(seconds=89100/10**9)*bcr_count
        self.bcr_count = bcr_count
        self.triggers = []
    
    def __str__(self):
        return f"""BCR({self.hitTime}*25/32 ns, {self.absTime}, {self.event_time},{self.bcr_count}) \n    Triggers: {[str(t) for t in self.triggers]})"""

    def add_trigger(self, trigger):
        trigger.absTime = self.absTime+(trigger.hitTime-self.hitTime)*datetime.timedelta(seconds=(25/32)/10**9)
        trigger.bcId = (trigger.hitTime-self.hitTime)/32.
        self.triggers.append(trigger)


class Trigger():
    def __init__(self, hitTime):
        self.hitTime = hitTime
        self.bcId = 0
        self.absTime = 0
    
    def __str__(self):
        return f"Trigger({self.hitTime}*25/32 ns, {self.absTime})"


class AtlasAnalyser():
    def __init__(self, file, trigger_channel, bcr_channel = 0):
        self.fReader = rawFileReader.fileReader(file)
        self.data = []
        self.trigger_channel = trigger_channel
        self.bcr_channel = bcr_channel
        
    def printBCR(self, amount_of_bcr):
        i = 0
        while i < amount_of_bcr:
            print(self.data[i])
            i += 1

    def printEvents(self, amount_of_events):
        last_event_time = datetime.datetime(1960,2,3) 
        evt_num = -1 #offset one
        bcr_num = 0
        while evt_num < amount_of_events: #i dont like the code :(
            if self.data[bcr_num].event_time != last_event_time:
                evt_num += 1
                last_event_time = self.data[bcr_num].event_time
            if evt_num < amount_of_events:
                print(self.data[bcr_num])
            if bcr_num+1 < len(self.data):
                print("Check:", self.data[bcr_num+1].hitTime-self.data[bcr_num].hitTime)
                bcr_num += 1
            else:
                evt_num = amount_of_events
    
    def getTDC5Data(self, amount_of_events):
        i = 0
        initial_time = None
        previous_event_time = None
        while i < amount_of_events:

            if not self.fReader.readBlock():
                print("Finished Reading File")
                return
            if(self.fReader.hasEvents()):
                tdc5Reads = self.fReader.getTDCFiveEvents()
                if not tdc5Reads:
                    continue
                else:
                    i += 1

                if i == 1:
                    previous_event_time = tdc5Reads[0][0]
                    continue

                bcr_count = 1
                tdcEvent = tdc5Reads[0] #actual pile of data (timestamp, data)
                current_bcr = BCR(-1, tdcEvent[0], bcr_count)
                if not initial_time:
                    initial_time = tdcEvent[0]
                    print("Initial Time:", initial_time)
                for word in tdcEvent[1]:
                    tdcChannel = (word>>24)&0x7f
                    tdcHitTime = word&0xfffff
                    if tdcChannel == self.bcr_channel:
                        if tdcHitTime != current_bcr.hitTime: #avoid repetition
                            self.data.append(current_bcr)
                            # looking at the event time before
                            current_bcr = BCR(tdcHitTime, previous_event_time, bcr_count)
                            bcr_count += 1
                    elif tdcChannel == self.trigger_channel:
                        if not current_bcr:
                            print("No BCR found for trigger")
                            continue
                        if not current_bcr.triggers:
                            current_bcr.triggers.append(Trigger(tdcHitTime))
                        else:
                            if tdcHitTime != current_bcr.triggers[-1].hitTime:
                                current_bcr.add_trigger(Trigger(tdcHitTime))
                    else:
                        pass
                        #print("Unknown Channel", tdcChannel)
                
        self.data = self.data[1:] #remove the first empty BCR
        return self.data

    def getAtlasData(self):
        infile = '/eos/user/m/mireveri/anubis/athena/run/ATLASData/data-ANALYSIS/datafiles.root'
        testFile = uproot.open(infile)
        anaTree = testFile["analysis"]
        feats = [branch.name for branch in anaTree.branches]
        evtArr = anaTree.arrays(feats,library="pd")
        evtArr
        thisEvt = 9580
        print(datetime.datetime.fromtimestamp(evtArr['TimeStamp'][thisEvt][0]),evtArr['TimeStampNS'][thisEvt][0],evtArr['BCID'][thisEvt][0])