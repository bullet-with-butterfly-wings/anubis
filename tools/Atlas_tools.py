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
from tqdm import tqdm
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
        self.timeStamp = datetime.datetime.timestamp(event_time) #+bcr_count*89100/10**9
        self.bcr_count = bcr_count #Do I needed??
        self.cycle = 0
        self.triggers = []
        self.error = error
    
    def __str__(self):
        return f"""BCR({self.hitTime}*25/32 ns, {readTimeStamp(self.timeStamp)}, {self.event_time},{self.bcr_count}, {self.error}) \n    Triggers: {[str(t) for t in self.triggers]})"""

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
        self.TDC5Reads = []
        self.anubis_data = []
        self.atlas_data = []
        
    def printBCR(self, a,b):
        i = a
        while i < b:
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
        last = 0
        current = 0
        bad_stage = False
        counter = 0
        period = 1_048_575
        interval = 114_048
        cycle = 0
        with tqdm(total=len(self.anubis_data), desc=f"Correcting", unit='BCR') as pbar:
            while current < len(self.anubis_data):
                now_bcr = self.anubis_data[current]
                
                if current != 0 and now_bcr.hitTime - self.anubis_data[current-1].hitTime < -2e5: #magic number
                    cycle += 1
                now_bcr.cycle = cycle

                if now_bcr.error and not bad_stage:
                    #print("Bad stage started:", current)
                    bad_stage = True
                    last = current - 1
                    triggers = [] #handle triggers later
                good_bcr = self.anubis_data[last]
                
                if bad_stage:
                    #print("In bad stage:", current)
                    #print(now_bcr)
                    delta = now_bcr.hitTime - good_bcr.hitTime
                    if delta < 0:
                        delta += period
                    #print("Delta:", delta)
                    if  delta % interval < 50 or interval - (delta % interval) < 50: #no overflow
                        prob_bcr_count  = round(delta / interval)
                        obs_bcr_count = current - last
                        for bad_bcr in range(last+1, current): #remove and replace
                            del self.anubis_data[last+1]
                        #print(self.anubis_data[last-4:current+4])
                        for new_bcr in range(1,prob_bcr_count):
                            self.anubis_data.insert(last+new_bcr, BCR((good_bcr.hitTime+(new_bcr)*interval) % period, good_bcr.event_time, good_bcr.bcr_count)) # will need to think about other data
                            for trig in triggers:
                                previous = (good_bcr.hitTime+(new_bcr)*interval) % period
                                next = (good_bcr.hitTime+(new_bcr+1)*interval) % period
                                overflow = False
                                if previous > next:
                                    next += period
                                    overflow = True
                                if last < trig.hitTime < next or (overflow and previous < trig.hitTime + period < next):
                                    self.anubis_data[last+new_bcr].add_trigger(trig)
                                    triggers.remove(trig)
                                else:
                                    break
                        #print("Rerranged")
                        #print(triggers)
                        #self.printBCR(last,last+prob_bcr_count)
                        #print("Rerranged")
                        now_bcr.error = False
                        bad_stage = False
                        #print(self.anubis_data[last-4:current+4])

                        #bcr_number += prob_bcr_count-obs_bcr_count+1
                        #print("Sorted:", current)
                        #print("Events added:", prob_bcr_count-obs_bcr_count)
                        current += prob_bcr_count-obs_bcr_count
                    else:
                        triggers += now_bcr.triggers
                        if current - last > 9:
                            #print("Give up")
                            bad_stage = False
                            self.anubis_data[current+1].error = False

                    if not now_bcr.error and bad_stage:
                        #print("Fixed")
                        bad_stage = False
                pbar.update(1)
                current += 1
        return counter
                   
                   

    def _readingRoutine(self, tdc5Reads,trigger_channel, bcr_channel, previous_event_time, previous_last_bcr):
        bcr_count = 0
        tdcEvent = tdc5Reads[0] #actual pile of data (timestamp, data)
        current_bcr = None
        counter = 0
        period = 1_048_575
        interval = 114_048
        for word in tdcEvent[1]:
            tdcChannel = (word>>24)&0x7f
            tdcHitTime = word&0xfffff
            if tdcChannel == bcr_channel:
                if not current_bcr: #avoid repetition
                    delta = tdcHitTime - previous_last_bcr
                else:
                    delta = tdcHitTime - current_bcr.hitTime
                if delta < 0:
                    delta += period

                if not current_bcr or current_bcr.hitTime != tdcHitTime: 
                    bcr_count += 1            
                    current_bcr = BCR(tdcHitTime, previous_event_time, bcr_count, error = abs(delta - interval) > 50)
                    self.anubis_data.append(current_bcr)
           
            elif tdcChannel == trigger_channel:
                if not current_bcr: #corresponds to the last BCR
                    self.anubis_data[-1].add_trigger(Trigger(tdcHitTime))
                    continue
                if not current_bcr.triggers:
                    current_bcr.add_trigger(Trigger(tdcHitTime))
                else:
                    if tdcHitTime != current_bcr.triggers[-1].hitTime: #no repetition
                        current_bcr.add_trigger(Trigger(tdcHitTime))
            else:
                pass
                #print("Unknown Channel", tdcChannel)
        a = current_bcr.hitTime if current_bcr else previous_last_bcr
        return  tdcEvent[0], a#current_bcr.hitTime
   
    def getTDC5Data(self,file, trigger_channel, bcr_channel = 0, amount_of_events=100, start = None, end = None):
        initial_time = None
        previous_event_time = None
        fromPKL = (file[-4:] == ".pkl")
        if fromPKL:
            with open(file, 'rb') as inp:
                tdc5Events = pickle.load(inp)
        else:
            if not self.fReader or self.anubis_file != file: 
                print("New")
                self.fReader = rawFileReader.fileReader(file)
                self.anubis_file = file
        if start and end:
            total = end-start
            units = "seconds"
        else:
            total = amount_of_events
            units = "Events"  
        
        i = 0         
        with tqdm(total=total, desc=f"Reading", unit=units) as pbar: 
            while i < total:
                if not fromPKL:
                    if not self.fReader.readBlock():
                        raise EOFError("You have reached the end of the file")
                    
                    if(self.fReader.hasEvents()):
                        tdc5Reads = self.fReader.getTDCFiveEvents()
                        if not tdc5Reads:
                            continue
                    else:
                        continue
                else:
                    tdc5Reads = tdc5Events[i]

                if i == 0:
                    previous_event_time = tdc5Reads[0][0]
                    previous_last_bcr = 620958 #the ambiguity
                    i += 1
                    continue
                i += 1
                previous_event_time, previous_last_bcr = self._readingRoutine(tdc5Reads, trigger_channel, bcr_channel, previous_event_time, previous_last_bcr)
                pbar.update(self.anubis_data[-1].timeStamp - pbar.n)
            #print(previous_event_time)
        return self.anubis_data

    def saveTDC5Data(self, file, storage, start = None, end = None):
        if not self.fReader or self.anubis_file != file: 
            print("New fReader")
            self.fReader = rawFileReader.fileReader(file)
            self.anubis_file = file

        initial_timestamp = 0
        while initial_timestamp == 0:
            if not self.fReader.readBlock():
                raise EOFError("You have reached the end of the file")                    
            if(self.fReader.hasEvents()):
                tdc5Reads = self.fReader.getTDCFiveEvents()
                if not tdc5Reads:
                    continue
            else:
                continue
            print(tdc5Reads[0][0])
            initial_timestamp = datetime.datetime.timestamp(tdc5Reads[0][0])
        print(initial_timestamp)
        current_timestamp = initial_timestamp
        with tqdm(total=round(start-initial_timestamp), desc=f"Skipping", unit="seconds") as pbar: 
            while current_timestamp < start:
                #self.fReader.skip_events(2)
                if not self.fReader.readBlock():
                    raise EOFError("You have reached the end of the file")
                if(self.fReader.hasEvents()):
                    tdc5Reads = self.fReader.getTDCFiveEvents()
                    #print(tdc5Reads)
                    if not tdc5Reads:
                        continue
                else:
                    continue
                current_timestamp = datetime.datetime.timestamp(tdc5Reads[0][0])
                #print(current_timestamp - pbar.n)
                pbar.update(current_timestamp - pbar.n)
               
        with tqdm(total=round(end - start), desc=f"Reading", unit="seconds") as pbar: 
            while current_timestamp < end:
                if not self.fReader.readBlock():
                    raise EOFError("You have reached the end of the file")
                
                if(self.fReader.hasEvents()):
                    tdc5Reads = self.fReader.getTDCFiveEvents()
                    if not tdc5Reads:
                        continue
                else:
                    continue
                self.TDC5Reads.append(tdc5Reads)
                current_timestamp = datetime.datetime.timestamp(tdc5Reads[0][0])
                pbar.update(round(current_timestamp - pbar.n))

        with open(storage, "wb") as f:
            pickle.dump(self.TDC5Reads, f)
            print("Saved")
        return self.TDC5Reads

    def getAtlasData(self, file):
        self.atlas_file = file
        testFile = uproot.open(self.atlas_file)
        anaTree = testFile["analysis"]
        feats = [branch.name for branch in anaTree.branches]
        evtArr = anaTree.arrays(feats,library="pd")
        self.atlas_data = evtArr

        self.atlas_data = self.atlas_data.sort_values(by=["TimeStamp", "TimeStampNS"])
        return self.atlas_data
    
    def binary_search(self, anubis_pointer, hit_time, time_window):
        left = anubis_pointer
        right = len(self.anubis_data)
        while left < right:
            mid = (left + right) // 2
            if self.anubis_data[mid].timeStamp < hit_time - time_window:
                left = mid + 1
            else:
                right = mid
        return left

    def pairBCRwithEvents(self):
        anubis_pointer = 0
        atlas_pointer = 0
        time_window = 90e-6
        self.matches = []
        with tqdm(total=len(self.atlas_data), desc=f"Matching", unit='Events') as pbar:        
            for atlas_pointer in range(len(self.atlas_data)):
                hit_time = self.atlas_data.iloc[atlas_pointer]["TimeStamp"]+self.atlas_data.iloc[atlas_pointer]["TimeStampNS"]*1e-9
                self.matches.append([self.atlas_data.iloc[atlas_pointer],[]])
                i = anubis_pointer 
                while i < len(self.anubis_data):
                    bcr = self.anubis_data[i]
                    if bcr.timeStamp < hit_time - time_window:
                        anubis_pointer = self.anubis_data.index(bcr)
                        i += 1
                    if bcr.timeStamp < hit_time - 10*time_window: #if it is too far away do binary search
                        anubis_pointer = self.binary_search(anubis_pointer, hit_time, time_window)
                        i = anubis_pointer
                    if abs(bcr.timeStamp - hit_time) < time_window:
                        for trigger in bcr.triggers:
                            #if abs(trigger.bcId - self.atlas_data.iloc[atlas_pointer]["BCID"]) < 20:
                            self.matches[atlas_pointer][1].append(trigger)
                        i += 1
                    if bcr.timeStamp > hit_time + time_window:
                        break
                pbar.update(1)
        return self.matches
    
def BCRHistogram(data, atlas = False, plot = True):
    if atlas:
        hist = data["BCID"]
    else:
        hist = []
        problems = 0
        for bcr in data:
            for trigger in bcr.triggers:
                if trigger.bcId > 3564:
                    problems += 1
                else:
                    hist.append(round(trigger.bcId))
        #strange number of problems ~ 8193
    bins = [i for i in range(0,3565)]
    counts, _ = np.histogram(hist, bins=bins, density=True)
    if plot:
        bins = bins[:-1]
        plt.step(bins, counts)
        
        plt.xlabel('Time since last BCR (ns)')
        plt.xlim(0, 3565)
        plt.ylabel('Frequency')
        plt.title(f'Histogram of BCID {"Atlas" if atlas else "Anubis"}')
        plt.show()
    return counts