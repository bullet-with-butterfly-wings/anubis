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
import pandas
import pickle

#sys.path.insert(1, 'C://Users//Peter//OneDrive - University of Cambridge//Desktop//summer2//Osiris Temp//processing//python')
sys.path.append(os.path.join(sys.path[0], 'Osiris Temp', 'processing', 'python'))
import rawFileReader
import Reconstruction_tools as RTools
import Timing_tools as TTools
importlib.reload(RTools)
importlib.reload(TTools)

#scp -o "ProxyJump jd2052@gw.hep.phy.cam.ac.uk" jd2052@pcls.hep.phy.cam.ac.uk:/r04/atlas/revering/data/24_07/24_07_30/proAnubis_240730_2017.raw C:\Users\jony\Downloads
def readTimeStamp(timestamp):
        return datetime.datetime.fromtimestamp(timestamp)


class BCR():
    def __init__(self, hitTime, event_time, bcr_count, error = False):
        self.hitTime = hitTime
        self.event_time = event_time
        self.timeStamp = datetime.datetime.timestamp(event_time)+bcr_count*89100/10**9
        self.bcr_count = bcr_count
        self.triggers = []
        self.error = error
    
    def __str__(self):
        return f"""BCR({self.hitTime}*25/32 ns, {readTimeStamp(self.timeStamp)}, {self.event_time},{self.bcr_count}) \n    Triggers: {[str(t) for t in self.triggers]})"""

    def add_trigger(self, trigger):
        delta = (trigger.hitTime-self.hitTime)
        if delta < 0:
            delta += 1_048_575
            
        trigger.timeStamp = self.timeStamp+delta*(25/32)/10**9
        #print("Trigger", trigger.hitTime)
        #print("Bcr", self.hitTime)
        trigger.bcId = delta/32.

        #print("BCID", trigger.bcId)
        self.triggers.append(trigger)


class Trigger():
    def __init__(self, hitTime):
        self.hitTime = hitTime
        self.bcId = 0
        self.timeStamp = 0
    
    def __str__(self):
        return f"Trigger({self.hitTime}*25/32 ns, {self.timeStamp}, {self.bcId})"


class AtlasAnalyser():
    def __init__(self, ):
        self.anubis_file = None
        self.atlas_file = None
        self.fReader = None
        self.tDiff = []
        self.anubis_data = []
        self.atlas_data = []
        
    def printBCR(self, amount_of_bcr):
        i = 0
        while i < amount_of_bcr:
            print(self.anubis_data[i])
            i += 1


    def printEvents(self, amount_of_events):
        last_event_time = datetime.datetime(1960,2,3) 
        evt_num = -1 #offset one
        bcr_num = 0
        while evt_num < amount_of_events: #i dont like the code :(
            if self.anubis_data[bcr_num].event_time != last_event_time:
                evt_num += 1
                last_event_time = self.anubis_data[bcr_num].event_time
            if evt_num < amount_of_events:
                print(self.anubis_data[bcr_num])
            if bcr_num+1 < len(self.anubis_data):
                print("Check:", self.anubis_data[bcr_num+1].hitTime-self.anubis_data[bcr_num].hitTime)
                bcr_num += 1
            else:
                evt_num = amount_of_events

    def correctionBCR(self):
        #if difference < 10 => space evenly
        #
        last = 0
        current = 0
        bad_stage = False
        shift = 0
        histogram = []
        counter = 0
        while current < len(self.anubis_data):
            if self.anubis_data[current].error and not bad_stage:
               last = current - 1
               bad_stage = True
            elif not self.anubis_data[current].error and bad_stage: # transtition to good
                    exp = ((current-last)*114_048) % 1_048_575
                    actual = (self.anubis_data[current].hitTime - self.anubis_data[last].hitTime) % 1_048_575
                    histogram.append(current-last)
                    if not ((actual-exp) % 114_048 < 5 or abs((actual-exp) % 114_048 - 114_048) < 5) and not 91900 < (actual-exp) % 114_048 < 91907:
                        counter += current-last
                        print("Jump, ", current-last)
                        print("Expecting", exp)
                        print("Actual", actual)
                        print("Difference:", actual-exp)
                        print("Mod Difference:", (actual-exp) % 114_048)    
                        print("Error time:", self.anubis_data[current].event_time)
                    bad_stage = False
            current += 1
        plt.xlabel("Length")
        plt.ylabel("Frequency")
        plt.title("Distribution of the flagged region")
        plt.grid()
        plt.hist(histogram, bins = [i+0.5 for i in range(1, 25)], histtype="step", label=f'Length of the flagged region')
        plt.show()
        return counter
                   
                   

    def _readingRoutine(self, tdc5Reads,trigger_channel, bcr_channel, previous_event_time, previous_last_bcr):
        bcr_count = 1
        tdcEvent = tdc5Reads[0] #actual pile of data (timestamp, data)
        current_bcr = BCR(-1, previous_event_time, bcr_count)
        for word in tdcEvent[1]:
            tdcChannel = (word>>24)&0x7f
            tdcHitTime = word&0xfffff
            if tdcChannel == bcr_channel:
                if tdcHitTime != current_bcr.hitTime: #avoid repetition
                    if not current_bcr.hitTime == -1:
                        self.anubis_data.append(current_bcr)
                        bcr_count += 1
                        delta = tdcHitTime - current_bcr.hitTime 
                        if delta < 0:
                            delta += 1_048_575
                        self.tDiff.append(delta)
                    else:
                        delta = tdcHitTime - previous_last_bcr
                        if delta < 0:
                            delta += 1_048_575
                        #self.tDiff.append(delta)
                    current_bcr = BCR(tdcHitTime, previous_event_time, bcr_count, error = abs(delta - 114_048) > 5)
            
            elif tdcChannel == trigger_channel:
                if not current_bcr:
                    print("No BCR found for trigger")
                    continue
                if not current_bcr.triggers:
                    current_bcr.add_trigger(Trigger(tdcHitTime))
                else:
                    if tdcHitTime != current_bcr.triggers[-1].hitTime: #no repetition
                        current_bcr.add_trigger(Trigger(tdcHitTime))
            else:
                pass
                #print("Unknown Channel", tdcChannel)
        return  tdcEvent[0], current_bcr.hitTime

    def getTDC5Data(self,file, trigger_channel, bcr_channel = 0, amount_of_events=100, fromPKL = False, ):
        i = 0
        initial_time = None
        previous_event_time = None
        if fromPKL:
            with open(file, 'rb') as inp:
                tdc5Events = pickle.load(inp)
        else:
            if not self.fReader or self.anubis_file != file: 
                print("New")
                self.fReader = rawFileReader.fileReader(file)
                self.anubis_file = file
               
        while i < amount_of_events:
            if not fromPKL:
                if not self.fReader.readBlock():
                    print("Finished Reading File")
                    return
                if(self.fReader.hasEvents()):
                    tdc5Reads = self.fReader.getTDCFiveEvents()
                    if not tdc5Reads:
                        continue
            else:
                tdc5Reads = tdc5Events[i]

            i += 1
            if i == 1:
                previous_event_time = tdc5Reads[0][0]
                previous_last_bcr = 0 #the ambiguity
                continue

            previous_event_time, previous_last_bcr = self._readingRoutine(tdc5Reads, trigger_channel, bcr_channel, previous_event_time, previous_last_bcr)
            #print(previous_event_time)
        return self.anubis_data

    def getAtlasData(self, file):
        self.atlas_file = file
        testFile = uproot.open(self.atlas_file)
        anaTree = testFile["analysis"]
        feats = [branch.name for branch in anaTree.branches]
        evtArr = anaTree.arrays(feats,library="pd")
        self.atlas_data = evtArr

        self.atlas_data = self.atlas_data.sort_values(by=["TimeStamp", "TimeStampNS"])
        return self.atlas_data
    

    def pairBcIdWithAtlasHit(self):
        anubis_pointer = 0 #index of already checked data
        atlas_pointer = 0
        # you will need to align them before the run
        triggers = list(np.concatenate(np.asarray([bcr.triggers for bcr in self.anubis_data], dtype="object")))
        triggers.sort(key = lambda x: x.timeStamp)
        triggers = triggers[330_000:]
        window = 20_000e-9 # in ns on both sides
        dev = 5 #in bcId
        potential = []
        new_potential = []
        while anubis_pointer < len(triggers) and atlas_pointer < len(self.atlas_data["TimeStamp"]):
            atlas_hit = self.atlas_data.iloc[atlas_pointer]["TimeStamp"]+self.atlas_data.iloc[atlas_pointer]["TimeStampNS"]*1e-9
            if atlas_hit - window < triggers[anubis_pointer].timeStamp < atlas_hit + window:
                print("In window")
                print("Trigger:", triggers[anubis_pointer])
                print("Atlas:", self.atlas_data.iloc[atlas_pointer])
                #print("Anubis:", triggers[anubis_pointer].bcID)
                #if abs(triggers[anubis_pointer].bcID-self.atlas_data.iloc[atlas_pointer]["BCID"]) < dev:
                #    print("Same BCID")
                anubis_pointer += 1
            elif atlas_hit + window < triggers[anubis_pointer].timeStamp:
                print("Atlas increase:", (atlas_hit, triggers[anubis_pointer].timeStamp))
                atlas_pointer += 1
                potential.append(new_potential)
                new_potential = []
            elif atlas_hit - window > triggers[anubis_pointer].timeStamp:
                print("Anubis increase:",(atlas_hit, triggers[anubis_pointer].timeStamp))
                anubis_pointer += 1
        return potential
                    
    def pairBcIdInConstruction(self):
        anubis_pointer = 0 #index of already checked data
        atlas_pointer = 0
        # you will need to align them before the run
        triggers = list(np.concatenate(np.asarray([bcr.triggers for bcr in self.anubis_data], dtype="object")))
        triggers.sort(key = lambda x: x.timeStamp)
        triggers = triggers[330_000:]
        matches = []
        current_min = 1
        while anubis_pointer < len(triggers) and atlas_pointer < len(self.atlas_data["TimeStamp"]):
            atlas_hit = self.atlas_data.iloc[atlas_pointer]["TimeStamp"]+self.atlas_data.iloc[atlas_pointer]["TimeStampNS"]*1e-9
            if 0 < atlas_hit - triggers[anubis_pointer].timeStamp:
                current_min = min(current_min, atlas_hit - triggers[anubis_pointer].timeStamp)
                anubis_pointer += 1
            elif atlas_hit - triggers[anubis_pointer].timeStamp < 0:
                #print("Overshoot")
                current_min = min(current_min, abs(atlas_hit - triggers[anubis_pointer].timeStamp))
                matches.append((self.atlas_data.iloc[atlas_pointer], triggers[anubis_pointer], current_min))
                current_min = 1
                atlas_pointer += 1 
        return matches
    

