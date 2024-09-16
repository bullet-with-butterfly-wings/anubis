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
import Reconstruction_tools as RTools
import Timing_tools as TTools
importlib.reload(RTools)
importlib.reload(TTools)
import Analysis_tools as ATools
from scipy.optimize import curve_fit

        
class Timing_Analyser():
    def __init__(self, event_chunk, processsed_event, diffHists = None,
                 scDiffs = None, normDiffs = None):
        self.event_chunk = event_chunk
        self.processed_event = processsed_event
        self.tot_mean = [[[13.0 for eta in range(32)] for phi in range(64)] for rpc in range(6)]
        self.tot_std = [[[1.0 for eta in range(32)] for phi in range(64)] for rpc in range(6)]
        if diffHists:
            self.diffHists = diffHists
        else:
            self.diffHists = {i:[[hi.Hist(hi.axis.Regular(bins=150, start=-50, stop=50, name=f"rpc{i}etPhiDiff")) for etchan in range(32)] for phchan in range(64)] for i in range(6)}
        
        if scDiffs == None: #please erase
            self.scDiffs = [[0 for etchan in range(32)] for phchan in range(64)]
            self.normDiffs = [[0 for etchan in range(32)] for phchan in range(64)]
        else:
            self.scDiffs = scDiffs
            self.normDiffs = normDiffs
        self.residEtaLatest = []
        self.residPhiLatest = []
        self.count = [[] for _ in range(7)]
        self.rpc_involvement = {i: [0] * 6 for i in range(7)}
        
        
        
    def update_event(self, event_chunk, processed_event):
        self.event_chunk = event_chunk
        self.processedEvents = processed_event
        
    def read_tot_diffs(self):
        for event in self.event_chunk:
            Hits = TTools.tdcEventToRPCData(event,activeTDCs=[0,1,2,3,4])
            for rpc in range(6):
                minEtaHit = ATools.rpcHit(-1,1000,True, self.processed_event, rpc)
                minPhiHit = ATools.rpcHit(-1,1000,False, self.processed_event, rpc)
                for hit in Hits[rpc]:
                    if 150 < hit.time < 350:
                        if hit.eta and hit.time<minEtaHit.time:
                            minEtaHit = hit
                        elif hit.time<minPhiHit.time and not hit.eta:
                            minPhiHit = hit
                if minEtaHit.channel>-0.5 and minPhiHit.channel>-0.5:
                    self.diffHists[rpc][minPhiHit.channel][minEtaHit.channel].fill(minEtaHit.time-minPhiHit.time)
    
    def calculate_corrections(self):
        def gaus(x,a,x0,sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))
        for rpc in range(6):
            for phi in range(64):
                for eta in range(32):
                    if sum([1 for value in self.diffHists[rpc][phi][eta].values() if value != 0]) > 3:
                        try:
                            ppar, pcov = curve_fit(gaus, self.diffHists[rpc][phi][eta].axes.centers[0], self.diffHists[rpc][phi][eta].values(), p0=[1,13,5])
                            self.tot_mean[rpc][phi][eta] = ppar[1]
                            self.tot_std[rpc][phi][eta] = abs(ppar[2])
                        except:
                            values = []
                            for idx, value in enumerate(self.diffHists[rpc][phi][eta].values()):
                                for i in range(int(value)):
                                    values.append(self.diffHists[rpc][phi][eta].axes.centers[0][idx])
                            self.tot_mean[rpc][phi][eta] = np.mean(values)
                            self.tot_std[rpc][phi][eta] = np.std(values)
                            
                            #print(self.diffHists[rpc][phi][eta].values())
    
        return self.tot_mean, self.tot_std


        
    #so focking unsuable, I will leave it for now
    def Calculate_Residual_and_plot_TDC_Time_Diffs(self, pdf_filename='plots.pdf', max_itr=1):
        badPhi = {0: [], 1: [], 2: [31], 3: [0], 4: [19], 5: [31]}
        badEta = {0: [29, 30, 31], 1: [16, 20], 2: [], 3: [20, 31], 4: [], 5: []}
        
        rpcNames = {0: "Triplet Low", 1: "Triplet Mid", 2: "Triplet Top", 3: "Singlet", 4: "Doublet Low", 5: "Doublet Top"}
        evtCount = 0

        with matplotlib.backends.backend_pdf.PdfPages(pdf_filename) as pdf:
            for rpc in [0, 1, 2, 3, 4, 5]:
                slope = 0.15204446322001586
                offSet = 16.101626655986728
                for ph in range(64):
                    for et in range(32):
                        if sum(self.diffHists[rpc][ph][et].counts()) > 0:
                            self.scDiffs[ph][et] = sum( #mean
                                [thisVal * self.diffHists[rpc][ph][et].axes.centers[0][idx] for idx, thisVal in
                                enumerate(self.diffHists[rpc][ph][et])]) / sum(self.diffHists[rpc][ph][et].counts())
                            self.normDiffs[ph][et] = self.scDiffs[ph][et] + slope * (ph - et) - offSet

                fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
                etachannels = [x - 0.5 for x in range(33)]
                phichannels = [x - 0.5 for x in range(65)]
                etaHist = (self.scDiffs, np.array(phichannels), np.array(etachannels))
                zrange = [-5, 30]
                thisHist = hep.hist2dplot(etaHist, norm=colors.Normalize(zrange[0], zrange[1]))
                thisHist.cbar.set_label('$\eta$ hit time - $\phi$ hit time, Average (ns)', rotation=270, y=0.3, labelpad=23)
                plt.ylim(31.5, -0.5)
                plt.ylabel("Eta Channel")
                plt.xlabel("Phi Channel")
                ax.set_title(rpcNames[rpc])
                x_points = [-0.5, 64.5]
                y_points_dotted = [7.5, 15.5, 23.5]
                for y_point in y_points_dotted:
                    plt.plot(x_points, [y_point, y_point], 'k', linestyle='dotted')
                y_points_dashed = [-0.5, 31.5]
                x_points_dashed = [7.5, 15.5, 23.5, 31.5, 39.5, 47.5, 55.5]
                for x_point in x_points_dashed:
                    plt.plot([x_point, x_point], y_points_dashed, 'k', linestyle='dashed')
                pdf.savefig(fig)  # Save the figure to the PDF
                plt.close(fig)

                fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
                normHist = (self.normDiffs, np.array(phichannels), np.array(etachannels))
                zrange = [-4, 4]
                thisHist = hep.hist2dplot(normHist, norm=colors.Normalize(zrange[0], zrange[1]))
                thisHist.cbar.set_label('$\eta$ hit time - $\phi$ hit time, Average (ns)', rotation=270, y=0.3, labelpad=23)
                plt.ylim(31.5, -0.5)
                plt.ylabel("Eta Channel")
                plt.xlabel("Phi Channel")
                for y_point in y_points_dotted:
                    plt.plot(x_points, [y_point, y_point], 'k', linestyle='dotted')
                for x_point in x_points_dashed:
                    plt.plot([x_point, x_point], y_points_dashed, 'k', linestyle='dashed')
                ax.set_title(rpcNames[rpc] + ", Propagation Time Corrected")
                pdf.savefig(fig)  # Save the figure to the PDF
                plt.close(fig)

                phTimes = [0 for chan in range(64)]
                etTimes = [0 for chan in range(32)]

                phTimesBuffer = [0 for chan in range(64)]
                etTimesBuffer = [0 for chan in range(32)]
                testSpeed = 0.173
                offSet = 13.46
                itr = 0
                for phchan in range(64):
                    for etchan in range(32):
                        if etchan not in badEta[rpc]:
                            phTimes[phchan] += self.scDiffs[phchan][etchan] + (-etchan * testSpeed + phchan * testSpeed - offSet)
                    phTimes[phchan] /= (32. - len(badEta[rpc]))

                itr += 1

                for phchan in range(64):
                    for etchan in range(32):
                        if phchan not in badPhi[rpc]:
                            etTimes[etchan] += (
                                        self.scDiffs[phchan][etchan] + (-etchan * testSpeed + phchan * testSpeed - offSet) - phTimes[
                                    phchan]) / (64. - len(badPhi[rpc]))

                itr += 1

                while itr <= max_itr:
                    phTimesBuffer = phTimes.copy()
                    etTimesBuffer = etTimes.copy()
                    phTimes = [0 for chan in range(64)]
                    etTimes = [0 for chan in range(32)]
                    for phchan in range(64):
                        for etchan in range(32):
                            if etchan not in badEta[rpc]:
                                phTimes[phchan] += (self.scDiffs[phchan][etchan] + (
                                            -etchan * testSpeed + phchan * testSpeed - offSet) + etTimesBuffer[etchan]) / (
                                                        32. - len(badEta[rpc]))
                    for phchan in range(64):
                        for etchan in range(32):
                            if phchan not in badPhi[rpc]:
                                etTimes[etchan] += (self.scDiffs[phchan][etchan] + (
                                            -etchan * testSpeed + phchan * testSpeed - offSet) - phTimesBuffer[phchan]) / (
                                                        64. - len(badPhi[rpc]))

                    itr += 1
                self.residEtaLatest.append([time for time in etTimes])
                self.residPhiLatest.append([time for time in phTimes])
        return self.residEtaLatest, self.residPhiLatest
    
    def check_eta_trigger(self):
        etaHits = [[] for rpc in range(6)]
        skip_event = False
        for idx, event in enumerate(self.event_chunk):
            for tdc in range(5):
                if event.tdcEvents[tdc].qual != 0:
                    skip_event = True
                    break 
            skip_event = True
            if skip_event:
                continue 
            for tdc in range(5):
                for word in event.tdcEvents[tdc].words:
                    rpc, thisHit = ATools.tdcChanToRPCHit(word,tdc, self.processedEvents + idx)
                    if thisHit.channel == [0]:
                        continue
                    if thisHit.eta:
                        etaHits[rpc].append(thisHit)
        
        event_counts = {}

        for rpc_index, rpc_hits in enumerate(etaHits):
            unique_events_in_rpc = set()

            for hit in rpc_hits:
                unique_events_in_rpc.add(hit.event_num)

            for event_num in unique_events_in_rpc:
                if event_num not in event_counts:
                    event_counts[event_num] = 0
                event_counts[event_num] += 1

        for event_num, count in event_counts.items():
            if count > 6:
                print("WTF")
            if count >= 6:
                self.count[6].append(event_num)
            else:
                self.count[count].append(event_num)

            for rpc_index, rpc_hits in enumerate(etaHits):
                if any(hit.event_num == event_num for hit in rpc_hits):
                    if count >= 6:
                        self.rpc_involvement[6][rpc_index] += 1
                    else:
                        self.rpc_involvement[count][rpc_index] += 1

        # Check for failed events (count < 4)
        failed_events = [(event_num, count) for event_num, count in event_counts.items() if count < 4]

        if failed_events:
            return False, failed_events
        else:
            return True, None
        
    def plot_rpc_histogram(self):
        plt.figure(figsize=(10, 6))
        index = range(6)
        rpc_count = [0] * 6
        for count in range(6):
            rpc_count = list(map(sum, zip(rpc_count,self.rpc_involvement[count])))
        plt.bar(list(range(6)), rpc_count)

        plt.xlabel('RPC Number')
        plt.ylabel('Proportion of Events')
        plt.title('RPC Involvement')
        #plt.xticks([x + bar_width * 3 for x in index], index) 
        plt.legend()
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        plt.show()

    def plot_rpc_involvement_histogram(self):
        plt.figure(figsize=(10, 6))
        bar_width = 0.1 
        index = range(6)

        for i, count in enumerate(range(7)):
            total = sum(self.rpc_involvement[count])

            if total == 0:
                normalized = [0] * len(self.rpc_involvement[count])
            else:
                normalized = [rpc_count for rpc_count in self.rpc_involvement[count]]

            plt.bar([x + i * bar_width for x in index], normalized, bar_width, label=f'Count {count}')

        plt.xlabel('RPC Number')
        plt.ylabel('Proportion of Events')
        plt.title('Normalized RPC Involvement in Event Counts')
        plt.xticks([x + bar_width * 3 for x in index], index) 
        plt.legend()
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        plt.show()
        
    def plot_and_fit_tof_corrections(self):
        excludedEta = [31]
        excludedPhi = [19]
        fitPhChan = [0.5+chan for chan in range(64)]
        fitEtChan = [0.5+chan for chan in range(32)]
        fitDiffs = [row[:] for row in self.scDiffs]
        distance_per_phi_channel = 2.7625 #cm
        distance_per_eta_channel = 2.9844 #cm

        for phchan in range(64):
            for etchan in range(32):
                if phchan in excludedPhi:
                    fitDiffs[phchan][etchan] = self.scDiffs[phchan + 1][etchan]
                if etchan in excludedEta:
                    fitDiffs[phchan].pop(etchan)

        for etChan in excludedEta:
            fitEtChan.pop(etChan)
            
        # fitPhDist = [chan * distance_per_phi_channel for chan in fitPhChan]
        # fitEtDist = [chan * distance_per_eta_channel for chan in fitEtChan]

        def plane(xy, a, b):
            x, y = xy
            return (a * (x - y) + b).ravel()

        x, y = np.meshgrid(fitEtChan, fitPhChan)

        popt, pcov = curve_fit(plane, (x, y), np.array(fitDiffs).ravel(), p0=[0.175, 13])
        print(popt[0], popt[1])

        # Calculate fitted data and residuals
        data_fitted = plane((x, y), *popt)
        residuals = np.array(fitDiffs).ravel() - data_fitted
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((np.array(fitDiffs).ravel() - np.mean(np.array(fitDiffs).ravel()))**2)
        r_squared = 1 - (ss_res / ss_tot)

        print(f"R² value: {r_squared}")

        fig, ax = plt.subplots(1, 1)
        heatmap = ax.imshow(np.array(fitDiffs).reshape(64, len(fitEtChan)), cmap=plt.cm.jet, origin='lower',
                            extent=(0, max(fitEtChan), 0, max(fitPhChan)), vmin=0, vmax=30)
        plt.colorbar(heatmap, ax=ax)
        ax.contour(x, y, data_fitted.reshape(64, len(fitEtChan)), 16, colors='w')
        plt.xlabel('Eta channel')
        plt.ylabel('Phi channel')
        plt.title('Heatmap and Fitted Plane Contours')
        plt.show()
        
    def plot_residual(self):
        rpcNames = {0:"Triplet Low",1: "Triplet Mid", 2:"Triplet Top", 3:"Singlet",4:"Doublet Low",5:"Doublet Top"}
        fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
        phichannels = [x-0.5 for x in range(65)]

        for idx, rpc in enumerate([0,1,2,3,4,5]):
            plotPhiResids = self.residPhiLatest[idx].copy()
            plotPhiResids.append(plotPhiResids[-1])
            plt.step(phichannels,plotPhiResids,linewidth=3,label=rpcNames[rpc],where='post')
        yrange = ax.get_ylim()
        #ax.text(40, 0.8*yrange[1], "Slope: "+str(round(popt[0],3))+" ns/channel, Offset: "+str(round(popt[1],2))+" ns", fontsize=14,
        #               verticalalignment='top')
        ax.set_xlabel('$\phi$ Channel')
        ax.set_ylabel('Time Residual (ns)')
        # ax.set_ylim([-6,6])
        ax.set_xlim([-0.5,63.5])
        plt.legend()
        plt.show()
        fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
        etchannels = [x-0.5 for x in range(33)]
        for idx, rpc in enumerate([0,1,2,3,4,5]):
            plotEtaResids = self.residEtaLatest[idx].copy()
            plotEtaResids.append(plotEtaResids[-1])
            plt.step(etchannels,plotEtaResids,linewidth=3,label=rpcNames[rpc],where='post')
        yrange = ax.get_ylim()
        #ax.text(40, 0.8*yrange[1], "Slope: "+str(round(popt[0],3))+" ns/channel, Offset: "+str(round(popt[1],2))+" ns", fontsize=14,
        #               verticalalignment='top')
        ax.set_xlabel('$\eta$ Channel')
        ax.set_ylabel('Time Residual (ns)')
        # ax.set_ylim([-6,6])
        ax.set_xlim([-0.5,31.5])
        plt.legend()
        plt.show()     
    

