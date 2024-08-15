import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
hep.style.use([hep.style.ATLAS])
import sys
import os
import matplotlib.colors as colors
import matplotlib.backends.backend_pdf
import hist as hi
import importlib
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
        self.absTime = 0
        self.bcr_count = bcr_count # negative because timestamp at the end
        self.triggers = []
    
    def __str__(self):
        return f"""BCR({self.hitTime}*25/32 ns, {self.absTime}, {self.bcr_count}) \n    Triggers: {[str(t) for t in self.triggers]})"""

    def add_trigger(self, trigger):
        #set_the absolute time for the trigger
        self.triggers.append(trigger)

class Trigger():
    def __init__(self, hitTime):
        self.hitTime = hitTime
        self.absTime = 0
    
    def __str__(self):
        return f"Trigger({self.hitTime}*25/32 ns, {self.absTime})"

class AtlasAnalyser():
    def __init__(self, file, trigger_channel, bcr_channel = 0):
        self.fReader = rawFileReader.fileReader(file)
        self.data = []
        self.trigger_channel = trigger_channel
        self.bcr_channel = bcr_channel
        
    def printEvents(self, amount_of_bcr):
        i = 0
        while i < amount_of_bcr:
            print(self.data[i])
            i += 1

            
    def getTDC5Data(self, amount_of_events):
        i = 0
        initial_time = None
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
                
                bcr_count = 0
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
                            current_bcr = BCR(tdcHitTime, tdcEvent[0], bcr_count)
                            bcr_count += 1
                    elif tdcChannel == self.trigger_channel:
                        if not current_bcr:
                            print("No BCR found for trigger")
                            continue
                        if not current_bcr.triggers:
                            current_bcr.triggers.append(Trigger(tdcHitTime))
                        else:
                            if tdcHitTime != current_bcr.triggers[-1].hitTime:
                                current_bcr.triggers.append(Trigger(tdcHitTime))
                    else:
                        pass
                        #print("Unknown Channel", tdcChannel)
        self.data = self.data[1:] #remove the first empty BCR
        return self.data
