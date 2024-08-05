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
import Analysis_tools as ATools
import Reconstruction_tools as RTools
import Timing_tools as TTools
importlib.reload(RTools)
importlib.reload(TTools)
rpcHit = ATools.rpcHit
from scipy.optimize import curve_fit


class rpcCoincidence():
    def __init__(self, event_num, time_bin, hits):
        self.event_num = event_num
        self.time_bin = time_bin
        self.hits = hits

    def __str__(self):
        return f"rpcCoincidence(event_num={self.event_num}, time_bin={self.time_bin}, hits={self.hits})"
    

class Reconstructor():
    def __init__(self, event_chunk, processsed_event, tolerance = None, coincidence_window = 15, tof_correction = True):
        self.event_chunk = event_chunk
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]
        self.window_size = coincidence_window
        self.tof_correction = tof_correction
        self.processedEvents = processsed_event
        self.tol = [i for i in range(20)] if tolerance is None else tolerance
        self.dT = []
        self.recon = []
        
        
        
        self.possible_reconstructions = [0 for rppc in range(6)]
        self.successful_reconstructions = [[0 for i in range(len(self.tol))] for rpc in range(6)]
        self.successful_reconstructed_coords = {0:[[0 for etchan in range(32)] for phchan in range(64)],
                1:[[0 for etchan in range(32)] for phchan in range(64)], 
                2:[[0 for etchan in range(32)] for phchan in range(64)], 
                3:[[0 for etchan in range(32)] for phchan in range(64)],
                4:[[0 for etchan in range(32)] for phchan in range(64)],
                5:[[0 for etchan in range(32)] for phchan in range(64)]}
        self.failed_reconstructed_coords = {0:[[0 for etchan in range(32)] for phchan in range(64)],
                1:[[0 for etchan in range(32)] for phchan in range(64)], 
                2:[[0 for etchan in range(32)] for phchan in range(64)], 
                3:[[0 for etchan in range(32)] for phchan in range(64)],
                4:[[0 for etchan in range(32)] for phchan in range(64)],
                5:[[0 for etchan in range(32)] for phchan in range(64)]}
        self.possible_reconstructions_coords = {0:[[0 for etchan in range(32)] for phchan in range(64)],
                1:[[0 for etchan in range(32)] for phchan in range(64)], 
                2:[[0 for etchan in range(32)] for phchan in range(64)], 
                3:[[0 for etchan in range(32)] for phchan in range(64)],
                4:[[0 for etchan in range(32)] for phchan in range(64)],
                5:[[0 for etchan in range(32)] for phchan in range(64)]}
        
        
        self.eta_histogram = np.zeros(len(np.arange(-90.5, 91.5, 1)) - 1)
        self.phi_histogram = np.zeros(len(np.arange(-90.5, 91.5, 1)) - 1)
        self.solid_theta_histogram = np.zeros(len(np.arange(-180.5, 181.5, 1)) - 1)
        self.solid_phi_histogram = np.zeros(len(np.arange(-180.5, 181.5, 1)) - 1)
        self.tdcstatus = [True for tdc in range(5)]

        

    def populate_hits(self):
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]
        skip_event = False
        for idx, event in enumerate(self.event_chunk):
            for tdc in range(5):
                if event.tdcEvents[tdc].qual != 0:
                    skip_event = True
                    break 

            if skip_event:
                continue 
            for tdc in range(5):
                for word in event.tdcEvents[tdc].words:
                    rpc, thisHit = ATools.tdcChanToRPCHit(word,tdc, self.processedEvents + idx)
                    if thisHit.channel == [0]:
                        continue
                    if thisHit.eta:
                        self.etaHits[rpc].append(thisHit)

                    else:
                        self.phiHits[rpc].append(thisHit)
                        
    
    def update_event(self, event_chunk, processed_event):
        self.event_chunk = event_chunk
        self.processedEvents = processed_event
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]
        
    
    def make_cluster(self):
        coincident_hits = RTools.FindCoincidentHits(self.etaHits, self.phiHits, self.window_size, tof_correction=self.tof_correction)
        clustered = RTools.cluster(coincident_hits)
        return clustered
    
    
    def reconstruct_and_extrapolate(self, dataset, chi2_region = [0, 100]):
        # Ensure RPC is a list, even if it's a single integer
        if self.tdcstatus[3] == True:
            for rpc in range(6):
                for i, data in enumerate(dataset):
                    if RTools.count_entries(data) < 100:
                        E_recon = RTools.reconstruct_timed_Chi2_ByRPC(data, 3, rpc)
                        if E_recon:
                            if len(E_recon[2]) >= 5:
                                if E_recon[4] > chi2_region[0] and E_recon[4] < chi2_region[1]:
                                    # self.chi2.append(E_recon[4])
                                    # self.event_of_interest.append(E_recon)
                                    # Adding this check to see if other 5 RPCs are in reconstructed event.
                                    # This is necessary to ensure the reconstructed path is accurate.

                                    muon_coords = RTools.does_muon_hit_RPC(E_recon[0], E_recon[1], rpc)
                                    if muon_coords:
                                        self.possible_reconstructions[rpc] += 1
                                        self.possible_reconstructions_coords[rpc][int(muon_coords[0] / 2.7625)][int(muon_coords[1] / 2.9844)] += 1
                                        for idx, t in enumerate(self.tol):
                                            check = RTools.does_RPC_detect_muon(muon_coords, E_recon[7], t)
                                            if check:
                                                self.successful_reconstructions[rpc][idx] += 1
                                                self.successful_reconstructed_coords[rpc][int(muon_coords[0] / 2.7625)][int(muon_coords[1] / 2.9844)] += 1
                                            else:
                                                self.failed_reconstructed_coords[rpc][int(muon_coords[0] / 2.7625)][int(muon_coords[1] / 2.9844)] += 1


    def reconstruct_and_findtof(self, dataset, rpc_comparisons):
        for i, data in enumerate(dataset):
            if RTools.count_entries(data) < 100:
                E_recon = RTools.reconstruct_timed_Chi2_ByRPC(data, 3, -1, rpc_indicies=rpc_comparisons)
                if E_recon:
                    if len(E_recon[2]) >= 6:
                        self.dT.append(E_recon[5])
                        self.recon.append(E_recon)
                        
    def extract_angles_phi_eta_timed_DZ_modified(self, filtered_events, max_length=None, exact_length=False):
        angles_eta = []
        angles_phi = []
        delta_times = []
        dZ = []
        chi2_values = []
        solid_angle = []
        
        bin_size = 1
        eta_phi_bin_edges = np.arange(-90.5, 91.5, bin_size)
        solid_phi_bin_edges = np.arange(-180.5, 181.5, bin_size)

        for i, filtered_event in enumerate(filtered_events):
            if RTools.count_entries(filtered_event) < 100:
                result = RTools.reconstruct_timed_Chi2_modified(filtered_event, 3, max_length=max_length, exact_length=exact_length)
            else:
                result = None

            if result is not None:
                if result[1][2] > 0:
                    delta_times.append(result[5])
                    chi2_values.append(result[4])
                    dZ.append(result[6])

                    v_parr_eta = np.array([0, result[1][1], result[1][2]])
                    theta_eta = np.arccos(np.dot(v_parr_eta, [0, 0, 1]) / np.linalg.norm(v_parr_eta))

                    if theta_eta > np.pi / 2:
                        theta_eta = np.pi - theta_eta
                    if v_parr_eta[1] > 0:
                        theta_eta *= -1

                    angles_eta.append(theta_eta)

                    v_parr_phi = np.array([result[1][0], 0, result[1][2]])
                    theta_phi = np.arccos(np.dot(v_parr_phi, [0, 0, 1]) / np.linalg.norm(v_parr_phi))

                    if theta_phi > np.pi / 2:
                        theta_phi = np.pi - theta_phi
                    if v_parr_phi[0] < 0:
                        theta_phi *= -1

                    angles_phi.append(theta_phi)

                    magnitude = np.linalg.norm(result[1])
                    theta = np.arccos(result[1][2] / -magnitude)  # Polar angle
                    phi = np.arctan2(result[1][1], result[1][0])  # Azimuthal angle

                    # Normalize phi to the range [-π, π]
                    if phi > np.pi:
                        phi -= 2 * np.pi
                    elif phi < -np.pi:
                        phi += 2 * np.pi

                    # Convert angles to degrees
                    theta_deg = np.degrees(theta)
                    
                    phi_deg = np.degrees(phi)

                    # Append the solid angle (theta, phi) to self.angles_solid
                    solid_angle.append((theta_deg, phi_deg))

        angles_eta_degrees = [x * (180 / np.pi) for x in angles_eta]
        angles_phi_degrees = [x * (180 / np.pi) for x in angles_phi]
        angles_solid_theta_degrees = [x for x, y in solid_angle]
        angles_solid_phi_degrees = [y for x, y in solid_angle]

        self.eta_histogram += np.histogram(angles_eta_degrees, bins=eta_phi_bin_edges)[0]
        self.phi_histogram += np.histogram(angles_phi_degrees, bins=eta_phi_bin_edges)[0]
        self.solid_theta_histogram += np.histogram(angles_solid_theta_degrees, bins=solid_phi_bin_edges)[0]
        self.solid_phi_histogram += np.histogram(angles_solid_phi_degrees, bins=solid_phi_bin_edges)[0]
        
    def apply_systematic_correction(self, residEta, residPhi):
        for rpc in range(6):
            for i, etahit in enumerate(self.etaHits[rpc]):
                self.etaHits[rpc][i].time += residEta[rpc][etahit.channel]
            for j , phihit in enumerate(self.phiHits[rpc]):
                self.phiHits[rpc][j].time += residPhi[rpc][phihit.channel]
                
                
    def plot_tof_offset(self, rpc_comparison):
        tof = [[] for _ in range(6)]
        for dT in self.dT:
            for i in range(len(rpc_comparison)):
                tof[i].append(dT[i])
        Test_coord = [[] for _ in range(6)]
        for recon in self.recon:
            distance_per_phi_channel = 2.7625 #cm
            distance_per_eta_channel = 2.9844 #cm
            for i in range(6):
                recon[2][i][0] = int(recon[2][i][0] / distance_per_phi_channel)
                recon[2][i][1] = int(recon[2][i][1] / distance_per_eta_channel)
                Test_coord[i].append([recon[2][i][0], recon[2][i][1]])
                
        Coordinate = {0:[[[] for etchan in range(32)] for phchan in range(64)],
            1:[[[] for etchan in range(32)] for phchan in range(64)],
            2:[[[] for etchan in range(32)] for phchan in range(64)],
            3:[[[] for etchan in range(32)] for phchan in range(64)], 
                    4:[[[] for etchan in range(32)] for phchan in range(64)],
                        5:[[[] for etchan in range(32)] for phchan in range(64)]}

        scDiffs = [[0 for etchan in range(32)] for phchan in range(64)]
        normDiffs = [[0 for etchan in range(32)] for phchan in range(64)]
        rpcNames = {0:"Triplet Low",1: "Triplet Mid", 2:"Triplet Top", 3:"Singlet",4:"Doublet Low",5:"Doublet Top"}
        for i in range(len(tof[0])):
            for j in range(5):
                Coordinate[j + 1][Test_coord[j + 1][i][0]][Test_coord[j + 1][i][1]].append(tof[j][i])

                
        for rpc in [1,2,3,4,5]:
            for ph in range(64):
                for et in range(32):
                    # if sum(width[rpc][ph][et].counts())>0:
                        scDiffs[ph][et]=np.mean(Coordinate[rpc][ph][et])
            
            fig, ax = plt.subplots(1, figsize=(16, 8), dpi=100)
            etachannels = [x-0.5 for x in range(33)]
            phichannels = [x-0.5 for x in range(65)]
            etaHist = (scDiffs,np.array(phichannels),np.array(etachannels))
            zrange = [-20,30]
            thisHist = hep.hist2dplot(etaHist,norm=colors.Normalize(zrange[0],zrange[1]))
            thisHist.cbar.set_label('tof with 0 as reference', rotation=270, y=0.1,labelpad=23)
            plt.ylim(31.5,-0.5)
            plt.ylabel("Eta Channel")
            plt.xlabel("Phi Channel")
            ax.set_title(rpcNames[rpc])
            x_points = [-0.5, 64.5]
            y_points = [7.5, 15.5, 23.5]
            for y_point in y_points:
                plt.plot(x_points, [y_point,y_point], 'k', linestyle='dotted')
            y_points = [-0.5, 31.5]
            x_points = [7.5,15.5, 23.5, 31.5, 39.5, 47.5, 55.5]
            for x_point in x_points:
                plt.plot([x_point,x_point], y_points, 'k', linestyle='dashed')
            plt.show()
    

        
        
class Timing_Analyser():
    def __init__(self, event_chunk, processsed_event, diffHists = None,
                 scDiffs = None, normDiffs = None):
        self.event_chunk = event_chunk
        self.processed_event = processsed_event
        self.totDiffs = {0:[[0 for etchan in range(32)] for phchan in range(64)],
                1:[[0 for etchan in range(32)] for phchan in range(64)], 
                2:[[0 for etchan in range(32)] for phchan in range(64)], 
                3:[[0 for etchan in range(32)] for phchan in range(64)],
                4:[[0 for etchan in range(32)] for phchan in range(64)],
                5:[[0 for etchan in range(32)] for phchan in range(64)]}
        self.nDiffs = {0:[[0 for etchan in range(32)] for phchan in range(64)],
                    1:[[0 for etchan in range(32)] for phchan in range(64)], 
                    2:[[0 for etchan in range(32)] for phchan in range(64)], 
                    3:[[0 for etchan in range(32)] for phchan in range(64)],
                    4:[[0 for etchan in range(32)] for phchan in range(64)],
                    5:[[0 for etchan in range(32)] for phchan in range(64)]}
        if diffHists:
            self.diffHists = diffHists
        else:
            self.diffHists = {0:[[hi.Hist(hi.axis.Regular(bins=376, start=-150.4, stop=150.4, name="rpc0etPhiDiff")) for etchan in range(32)] for phchan in range(64)],
                        1:[[hi.Hist(hi.axis.Regular(bins=376, start=-150.4, stop=150.4, name="rpc1etPhiDiff")) for etchan in range(32)] for phchan in range(64)],
                        2:[[hi.Hist(hi.axis.Regular(bins=376, start=-150.4, stop=150.4, name="rpc2etPhiDiff")) for etchan in range(32)] for phchan in range(64)],
                        3:[[hi.Hist(hi.axis.Regular(bins=376, start=-150.4, stop=150.4, name="rpc3etPhiDiff")) for etchan in range(32)] for phchan in range(64)],
                        4:[[hi.Hist(hi.axis.Regular(bins=376, start=-150.4, stop=150.4, name="rpc4etPhiDiff")) for etchan in range(32)] for phchan in range(64)],
                        5:[[hi.Hist(hi.axis.Regular(bins=376, start=-150.4, stop=150.4, name="rpc5etPhiDiff")) for etchan in range(32)] for phchan in range(64)]}
        if scDiffs == None:
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
        
    def readTDCTimeDiffs(self):
        for event in self.event_chunk:
            Hits = TTools.tdcEventToRPCData(event,activeTDCs=[0,1,2,3,4])
            for rpc in [0,1,2,3,4,5]:
                minEtaHit = rpcHit(-1,1000,True, self.processed_event, rpc)
                minPhiHit = rpcHit(-1,1000,False, self.processed_event, rpc)
                for hit in Hits[rpc]:
                    if hit.time>150 and hit.time<300:
                        if hit.eta and hit.time<minEtaHit.time:
                            minEtaHit = hit
                        elif hit.time<minPhiHit.time and not hit.eta:
                            minPhiHit = hit
                if minEtaHit.channel>-0.5 and minPhiHit.channel>-0.5:
                    self.totDiffs[rpc][minPhiHit.channel][minEtaHit.channel] = self.totDiffs[rpc][minPhiHit.channel][minEtaHit.channel]+minEtaHit.time-minPhiHit.time
                    self.nDiffs[rpc][minPhiHit.channel][minEtaHit.channel] = self.nDiffs[rpc][minPhiHit.channel][minEtaHit.channel]+1
                    self.diffHists[rpc][minPhiHit.channel][minEtaHit.channel].fill(minEtaHit.time-minPhiHit.time)

    
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
                            self.scDiffs[ph][et] = sum(
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
                normalized = [rpc_count / total for rpc_count in self.rpc_involvement[count]]

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
    

