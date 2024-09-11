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
import gc
importlib.reload(RTools)
importlib.reload(TTools)

#scp -o "ProxyJump jd2052@gw.hep.phy.cam.ac.uk" jd2052@pcls.hep.phy.cam.ac.uk:/r04/atlas/revering/data/24_07/24_07_30/proAnubis_240730_2017.raw C:\Users\jony\Downloads
def readTimeStamp(timestamp):
        return datetime.datetime.fromtimestamp(timestamp)


class BCR():
    def __init__(self, hitTime, event_time, error = False):
        self.hitTime = hitTime
        self.event_time = event_time
        self.timeStamp = datetime.datetime.timestamp(event_time) #+bcr_count*89100/10**9
        self.triggers = []
        self.error = error
    
    def __str__(self):
        return f"""BCR({self.hitTime}*25/32 ns, {readTimeStamp(self.timeStamp)}, {self.event_time}, {self.error}) \n    Triggers: {[str(t) for t in self.triggers]})"""

    def add_trigger(self, trigger):
        delta = (trigger.hitTime-self.hitTime)
        if delta < 0:
            delta += 1_048_575

        if delta/32 > 3564 and not self.error:
            pass
        trigger.timeStamp = self.timeStamp+delta*(25/32)/10**9
        #print("Trigger", trigger.hitTime)
        #print("Bcr", self.hitTime)
        trigger.bcId = delta/32.
        """
        if trigger.bcId > 3564:
            print("Problem")
            print("Trigger", trigger.hitTime)
            print("Bcr", self.hitTime)
            print(self.error)
            print(self.triggers)
        #print("BCID", trigger.bcId)
        """
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
        
    def print_bcr(self, a,b):
        i = a
        while i < b:
            print(self.anubis_data[i])
            i += 1


    def print_events(self, amount_of_events):
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
    
    def correction_bcr(self):
        buffer = self.anubis_data.copy()
        interval = 100_000
        with tqdm(total=len(buffer), desc=f"Spliting", unit='BCR') as pbar:
            while buffer:
                storage = f"data/bcr/bcr{pbar.n}.pkl"
                with open(storage, "wb") as f:
                    pickle.dump(buffer[:interval], f)
                buffer = buffer[interval:]
                pbar.update(interval)
        
        data_list = sorted([f for f in os.listdir("data//bcr") if os.path.isfile(os.path.join("data//bcr", f))], key=lambda x: int(x[3:-4]))
        with tqdm(total=len(data_list), desc=f"Correcting", unit='files') as pbar:
            for data in data_list:
                with open(f"data//bcr//{data}", "rb") as f:
                    bcrs = pickle.load(f)
                bcrs = self.correction_bcr_routine(bcrs)
                with open(f"data//bcr//{data}", "wb") as f:
                    pickle.dump(bcrs, f)
                pbar.update(1)

        for data in data_list: #necessary because of garbage collector ig
            with open(f"data//bcr//{data}", "rb") as f:
                self.anubis_data += pickle.load(f)
        
        with open(f"data//corrected_bcr.pkl", "wb") as f:
            pickle.dump(self.anubis_data, f)
        print("Error rate:", sum([bcr.error*1 for bcr in self.anubis_data])/len(self.anubis_data))
        return self.anubis_data




    def correction_bcr_routine(self, bcrs):
        last = 0
        current = 0
        bad_stage = False
        period = 1_048_575
        interval = 114_048
        while current < len(bcrs):
            now_bcr = bcrs[current]
            if now_bcr.error and not bad_stage:
                bad_stage = True
                last = current - 1
                good_bcr = bcrs[last]
                triggers = good_bcr.triggers.copy() #handle triggers later
                good_bcr.triggers = []
            
            if bad_stage:
                delta = now_bcr.hitTime - good_bcr.hitTime
                if delta < 0:
                    delta += period
                if  delta % interval < 50 or interval - (delta % interval) < 50: #no overflow
                    prob_bcr_count  = round(delta / interval)
                    obs_bcr_count = current - last
                    for bad_bcr in range(last+1, current): #remove and replace
                        del bcrs[last+1]
                    for new_bcr in range(0,prob_bcr_count):
                        if new_bcr != 0:
                            bcrs.insert(last+new_bcr, BCR((good_bcr.hitTime+(new_bcr)*interval) % period, good_bcr.event_time))
                        for trig in triggers:
                            previous = (good_bcr.hitTime+(new_bcr)*interval) % period
                            next = (good_bcr.hitTime+(new_bcr+1)*interval) % period
                            overflow = False
                            if previous > next:
                                next += period
                                overflow = True
                            if previous < trig.hitTime < next or (overflow and previous < trig.hitTime + period < next):
                                bcrs[last+new_bcr].add_trigger(trig)
                                triggers.remove(trig)
                            else:
                                break
                    now_bcr.error = False
                    bad_stage = False
                    current += prob_bcr_count-obs_bcr_count
                else:
                    triggers += now_bcr.triggers
                    if current - last > 9:
                        bad_stage = False
                        bcrs[current+1].error = False

                if not now_bcr.error and bad_stage:
                    bad_stage = False
            current += 1
        return bcrs
                   
                   

    def _reading_routine(self, tdc5Reads,trigger_channel, bcr_channel, previous_event_time, previous_last_bcr):
        tdcEvent = tdc5Reads #actual pile of data (timestamp, data)
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
                    current_bcr = BCR(tdcHitTime, previous_event_time, error = abs(delta - interval) > 50)
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
   
    def get_tdc5_data(self,file, trigger_channel, bcr_channel = 0, amount_of_events=100, from_bcr = False):
        initial_time = None
        previous_event_time = None
        fromPKL = (file[-4:] == ".pkl")
        if from_bcr: #everything is already in the file
            with open(file, 'rb') as inp:
                self.anubis_data = pickle.load(inp)
            return self.anubis_data[:amount_of_events]
        if fromPKL:
            with open(file, 'rb') as inp:
                tdc5Events = pickle.load(inp)
        else:
            if not self.fReader or self.anubis_file != file: 
                print("New")
                self.fReader = rawFileReader.fileReader(file)
                self.anubis_file = file
    
        with tqdm(total=amount_of_events, desc=f"Reading", unit="Events") as pbar: 
            i = 0     
            while i < amount_of_events:
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
                    initial_time = tdc5Reads[0]
                    previous_event_time = tdc5Reads[0]
                    previous_last_bcr = 620958 #the ambiguity
                    i += 1
                    continue
                i += 1
                previous_event_time, previous_last_bcr = self._reading_routine(tdc5Reads, trigger_channel, bcr_channel, previous_event_time, previous_last_bcr)
                pbar.update(1)
            #print(previous_event_time)
        return self.anubis_data


    def get_atlas_data(self, file):
        self.atlas_file = file
        testFile = uproot.open(self.atlas_file)
        anaTree = testFile["analysis"]
        feats = [branch.name for branch in anaTree.branches]
        evtArr = anaTree.arrays(feats,library="pd")
        self.atlas_data = evtArr

        self.atlas_data = self.atlas_data.sort_values(by=["TimeStamp", "TimeStampNS"])
        return self.atlas_data
    
    #check
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

    def match_bcrs(self):
        anubis_pointer = 0
        atlas_pointer = 0
        time_window = 5e-3
        self.matches = []
        with tqdm(total=len(self.atlas_data), desc=f"Matching", unit='Events') as pbar:        
            for atlas_pointer in range(len(self.atlas_data)):
                hit_time = self.atlas_data.iloc[atlas_pointer]["TimeStamp"]+self.atlas_data.iloc[atlas_pointer]["TimeStampNS"]*1e-9
                self.matches.append([self.atlas_data.iloc[atlas_pointer].to_dict(), []])
                i = anubis_pointer 
                if atlas_pointer >= 34705:
                    self.print_bcr(i-4, i+4)
                    print("Atlas",hit_time)
                    print(self.atlas_data.iloc[atlas_pointer-1]["TimeStamp"]+self.atlas_data.iloc[atlas_pointer-1]["TimeStampNS"]*1e-9)
                while i < len(self.anubis_data):
                    bcr = self.anubis_data[i]
                    if atlas_pointer >= 34705:
                        print("New")
                        print(i)
                        print(bcr.timeStamp)
                        print(hit_time)
                        print("In:",abs(bcr.timeStamp - hit_time) < time_window)
                        print("High:",bcr.timeStamp > hit_time + time_window)
                        print("Low:", bcr.timeStamp < hit_time - time_window)
                        return self.matches
                    if bcr.timeStamp < hit_time - time_window:
                        """
                        anubis_pointer = self.anubis_data.index(bcr)
                        i += 1
                    if bcr.timeStamp < hit_time - 10*time_window: #if it is too far away do binary search
                        """
                        anubis_pointer = self.binary_search(anubis_pointer, hit_time, time_window)
                        i = anubis_pointer          
                        if atlas_pointer >= 34705:
                            print(i)
                        bcr = self.anubis_data[i-1]
                        if abs(bcr.timeStamp - hit_time) < time_window and not bcr.error:
                            for trigger in bcr.triggers:
                                if atlas_pointer >= 34705:
                                    print(trigger)
                                self.matches[atlas_pointer][-1].append(trigger)
                        if atlas_pointer >= 34705:
                            print("Done")
                    elif abs(bcr.timeStamp - hit_time) < time_window:
                        if not bcr.error:
                            for trigger in bcr.triggers:                            #if abs(trigger.bcId - self.atlas_data.iloc[atlas_pointer]["BCID"]) < 20:
                                self.matches[atlas_pointer][-1].append(trigger)
                        i += 1
                    elif bcr.timeStamp > hit_time + time_window:
                        if atlas_pointer >= 34705:
                            print("Break")
                        break
                pbar.update(1)
        return self.matches
    
    
def bcr_histogram(data, atlas = False, plot = True):
    if atlas:
        hist = data["BCID"]
    else:
        hist = []
        problems = 0
        for bcr in data:
            for trigger in bcr.triggers:
                if round(trigger.bcId) > 3564:
                    problems += 1
                else:
                    hist.append(round(trigger.bcId))
        #strange number of problems ~ 8193
    bins = [i-0.5 for i in range(0,3565)] #-0.5,0.5,...3653.5
    counts, _ = np.histogram(hist, bins=bins, density=True)
    if plot:
        bins = [i for i in range(0,3564)]
        plt.plot(bins, counts, color="orange")
        
        plt.xlabel('Time since last BCR (ns)')
        plt.xlim(0, 3565)
        plt.ylabel('Frequency')
        plt.title(f'Histogram of BCID {"Atlas" if atlas else "Anubis"}')
        plt.show()
    return counts

def positionn_filter_atlas(data, eta_func, phi_func):
    with tqdm(total=len(data.values)) as pbar:
        good_indices = []
        for i in range(len(data.values)):
            hits = data.iloc[i]
            for p in range(len(hits["mu_eta"])): #0.93, 0.83, 1.57
                if  eta_func(hits["mu_eta"][p]) and phi_func(hits["mu_phi"][p]):
                    good_indices.append(i)
                    break
            pbar.update(1)
    data = data.iloc[good_indices]
    return data
def beam_luminosity(data, plot = True):
    #beam luminosity
    times = data["TimeStamp"]
    hist, bins = np.histogram(times, bins=100)
    plt.plot(bins[:-1], hist)
    plt.show()

def eta_distribution(data, plot = True): #returns a list of etas
    hist_eta = []
    for i in range(len(data.values)):
        hits = data.iloc[i]
        for p in range(len(hits["mu_eta"])):
                hist_eta.append(hits["mu_eta"][p])

    if plot:
        plt.title("Density of muons by eta")
        plt.xlabel("Eta")
        plt.xlim(-3, 3)
        plt.hist(hist_eta, bins=100, density=True)
        plt.show()
    return hist_eta

"""
# set the times
analyser.atlas_data = analyser.atlas_data[analyser.atlas_data["TimeStamp"] > analyser.anubis_data[0].timeStamp]
print("Anubis", (round(analyser.anubis_data[0].timeStamp), round(analyser.anubis_data[-1].timeStamp)))
print("ATLAS",(analyser.atlas_data.iloc[0]["TimeStamp"], analyser.atlas_data.iloc[-1]["TimeStamp"]))
"""

def convert_matches(matches, best_offset):
    eta = []
    phi = []
    times = []
    heatmaps = [[0 for x in range(3564)] for y in range(3564)]
    matches = [x for x in matches if x[1]]
    for i in range(len(matches)):
        atlas_hits = matches[i][0]
        mu_eta = list(atlas_hits["mu_eta"])
        mu_phi = list(atlas_hits["mu_phi"])

        for p in range(len(mu_phi)): 
            #if 0.95 > mu_eta[p] > 0.8 and abs(mu_phi[p]) < 0.1:        
            if matches[i][1]:
                for trigger in matches[i][1]:
                    #if round(trigger.bcId+best_offset-matches[i][0]["BCID"]) < 3500:
                    #if abs(round(trigger.bcId-matches[i][0]["BCID"])) > 1000:
                    #    continue
                    if trigger.bcId > 3564:
                        continue
                    anubis_hit = round(trigger.bcId+best_offset) % 3564
                    
                    heatmaps[anubis_hit][matches[i][0]["BCID"]] += 1
                    #diffs.append(cdf[round(trigger.bcId+best_offset)] - cdf[matches[i][0]["BCID"]])
                    times.append(anubis_hit-matches[i][0]["BCID"])
                    if abs(anubis_hit-matches[i][0]["BCID"]) < 1:
                            eta.extend(mu_eta)
                            phi.extend(mu_phi)        

                    """
                    if round((trigger.bcId+best_offset) - matches[i][0]["BCID"]) < 0: #pod 
                        if abs(round((trigger.bcId+best_offset) - matches[i][0]["BCID"])) < 3564 - abs(round((trigger.bcId+best_offset) - matches[i][0]["BCID"])):
                            times.append(round((trigger.bcId+best_offset) - matches[i][0]["BCID"]))
                        else:
                            times.append(3564 + round((trigger.bcId+best_offset) - matches[i][0]["BCID"]))
                    else: #trigger nad tebou
                        if abs(round((trigger.bcId+best_offset) - matches[i][0]["BCID"])) < 3564 - abs(round((trigger.bcId+best_offset) - matches[i][0]["BCID"])):
                            times.append(round((trigger.bcId+best_offset) - matches[i][0]["BCID"]))
                        else:
                            times.append(-3564 + round((trigger.bcId+best_offset) - matches[i][0]["BCID"]))
                    """
                break
    return times, eta, phi, heatmaps

"""
# Define the step sizes for x and y
phi_step = 0.1  # Step size for x-axis
eta_step = 0.2  # Step size for y-axis

# Generate the ranges based on step size
eta_range = np.arange(-2.5, 2.5 + eta_step, eta_step)  # Range from -2.5 to 2.5
phi_range = np.arange(-np.pi, np.pi + phi_step, phi_step)  # Range from -pi to pi

# Generate the grid using meshgrid
phi_coordinates, eta_coordinates = np.meshgrid(phi_range, eta_range)
#print(phi_coordinates) #changes phi
#print(eta_coordinates) #changes eta

# Iterate over every single point in the grid
buffer = matches.copy()
objects = [[[] for i in range(phi_coordinates.shape[1])] for j in range(phi_coordinates.shape[0])]
with tqdm(total=len(buffer)) as pbar:
    for hit in matches:
        found = False
        pbar.update(1)
        for i in range(phi_coordinates.shape[0]):
            for j in range(phi_coordinates.shape[1]):
                phi_val = phi_coordinates[i, j]
                eta_val = eta_coordinates[i, j]
                #print(f"Point ({phi_val}, {eta_val})")
                for p in range(len(hit[0]["mu_eta"])):
                    if eta_val+eta_step/2 > hit[0]["mu_eta"][p] > eta_val-eta_step/2 and abs(hit[0]["mu_phi"][p] - phi_val) < phi_step/2:
                        for trigger in hit[1]:
                            if abs(round(trigger.bcId+best_offset-matches[i][0]["BCID"])) < 1000:
                                objects[i][j].append(round(trigger.bcId+best_offset-matches[i][0]["BCID"])) #add times
                        found = True
                        break
                if found:
                    break
            if found:
                break        

from collections import Counter
std_grid = 1000*np.ones((phi_coordinates.shape[0], phi_coordinates.shape[1]))

minimum = 10_000
# Iterate through the grid cells
for i in range(len(objects)):
    for j in range(len(objects[i])):
        if len(objects[i][j]) > 0:
            # Calculate the standard deviation of the numbers in this grid cell
            std_grid[i][j] = np.std(objects[i][j])#(Counter(objects[i][j]).most_common(1)[0][1]) if len(objects[i][j]) > 1 else 0
            if minimum > std_grid[i][j]:
                minimum = std_grid[i][j]
                print((i,j))
                print(phi_coordinates[i, j], eta_coordinates[i, j])
# Output the standard deviation array
print("Standard deviation for each grid cell:")
print(std_grid)

import matplotlib.pyplot as plt
import numpy as np

# Assuming std_grid is already calculated and phi_coordinates, eta_coordinates are defined

# Plotting a heatmap of std_grid
plt.imshow(std_grid, cmap='plasma', extent=[eta_coordinates.min(), eta_coordinates.max(),
                                            phi_coordinates.min(), phi_coordinates.max()],
           origin='lower', aspect='auto')

plt.colorbar(label='Standard Deviation')
plt.xlabel("Eta")
plt.ylabel("Phi")
plt.title("Heatmap of Standard Deviations")
plt.show()

    
"""
def generate_random_triggers(lenthg):
    with open("data/bcr_histogram.pkl", "rb") as f:
        density = pickle.load(f)
    found = 0
    random_bcr = []
    bound = max(density)
    while found < lenthg:
        bx = np.random.randint(0, 3564)
        value = np.random.uniform(0, bound)
        if value < density[bx]:
            random_bcr.append(bx)
            found += 1
    return random_bcr

def plot_offset(analyser):
    #Calculate the BCR offset
    anubis = bcr_histogram(analyser.anubis_data, atlas=False, plot = False)
    atlas = bcr_histogram(analyser.atlas_data, atlas=True, plot = False)

    window_size = 100
    #anubis = np.convolve(np.append(anubis, anubis[:window_size]), np.ones(window_size)/window_size, mode='valid')
    #atlas = np.convolve(np.append(atlas, atlas[:window_size]), np.ones(window_size)/window_size, mode='valid')
    error = 100
    errors = []

    for offset in range(-200,200):
        errors.append(np.sum((atlas - np.roll(anubis, offset))**2))
        if np.sum((atlas - np.roll(anubis, offset))**2) < error:
            error = np.sum((atlas - np.roll(anubis, offset))**2)
            best_offset = offset
    
    print("Best offset", best_offset)
    plt.plot(atlas, label="ATLAS")
    plt.plot(anubis, label="proANUBIS") #np.roll(anubis, best_offset)
    #plt.plot(anubis, label="Shifted proANUBIS") 
    plt.xlim(0, 3565)
    plt.ylabel("Counts (normalized)")
    plt.xlabel("BX")
    plt.title("proANUBIS and ATLAS Bunch Crossing Histograms - original")
    plt.legend()
    plt.show()

    plt.plot(list(range(-200,200)), errors, label="Errors")
    plt.xlabel("Offset")
    plt.ylabel("Sum of Square Errors")
    plt.title("Sum of Square Errors vs Offset")
    plt.vlines(best_offset,0, max(errors), color = "red", label=f"Best Offset: {best_offset}")
    plt.ylim(0,max(errors))
    plt.legend(loc = "lower left") #position

    plt.show()

def get_cdf(anubis_data):
    import scipy
    bcr_histogram = bcr_histogram(anubis_data)
    cdf = []
    total = 0
    for point in bcr_histogram:
        total += point
        cdf.append(total)

    print("CDF done")
    print(cdf)

    plt.title("CDF of proANUBIS BCR Histogram")
    plt.ylabel("CDF")
    plt.xlabel("BX")
    plt.plot(cdf)
    plt.show()
