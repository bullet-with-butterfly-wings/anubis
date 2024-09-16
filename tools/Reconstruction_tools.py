import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg', 'GTK3Agg', etc.
import mplhep as hep
hep.style.use([hep.style.ATLAS])
import sys
# import ANUBIS_triggered_functions as ANT
import matplotlib.backends.backend_pdf
from itertools import product
import numpy as np
# from scipy.stats import normpip install pillow
sys.path.insert(1, 'Osiris Temp\processing\python')
import Analysis_tools as ATools
from itertools import groupby
from itertools import product
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

time_range = (150,350)
RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.
tot_speed_factor = 0.15204446322001586 #(25/32) ns / channel

#time of flight, angles, 
class Track():
    def __init__(self, coordinates, uncertainties):
        self.coordinates = coordinates #t,x,y,z
        self.uncertainties = uncertainties #var_t, var_x, var_y
        self.centroid = np.mean(coordinates, axis=0)
        self.direction = None #proportional to velocity
        self.speed = 1 #c
        self.is_cosmic = False
        self.chi2 = 100 #arbitrary high value

    def __bool__(self):
        return True
    
    def fit(self):
        _, _, V = np.linalg.svd(self.coordinates-self.centroid)
        self.direction = V[0, :]
        """
        if self.direction[0] != 0:
            uncertainty = (self.uncertainties[0][1]*self.direction[1]**2 + self.uncertainties[0][2]*self.direction[2]**2)/np.linalg.norm(self.direction[1:])**4
            uncertainty += self.uncertainties[0][0]/self.direction[0]**2
            uncertainty = np.sqrt(uncertainty)
            
        """
        if (self.coordinates[-1][0] - self.coordinates[0][0]) != 0:
            self.speed = np.linalg.norm(self.coordinates[-1][1:]-self.coordinates[0][1:])/(self.coordinates[-1][0] - self.coordinates[0][0])
            self.speed *= (32/25*10**7)/2.99792458e8
        
        if self.direction[3] != 0:
            #self.speed = np.linalg.norm(self.direction[1:])/self.direction[0]*(32/25*10**7)
            #self.speed = (32/25*10**7)*np.linalg.norm(self.coordinates[-1][1:]-self.coordinate[0][1:])/(self.coordinates[-1][0] - self.coordinates[0][0])
            self.is_cosmic = self.direction[0]/self.direction[3] < 0
        self.calculate_chi2()   

    def calculate_chi2(self):
        chi2 = 0
        i = 0 #degrees of freeedom
        for idx, point in enumerate(self.coordinates):            
            i+=3
            if self.direction[3] == 0:
                chi2 = 100
                return chi2
            
            t = (point[3]-self.centroid[3])/self.direction[3] #parameter of the line
            for degree in range(1,3):
                ideal = self.centroid[degree] + t*self.direction[degree]
                real = point[degree]
                chi2+= (real-ideal)**2/self.uncertainties[idx][degree] 
        #check degrees of freedom
        doF = i - 6
        self.chi2 = chi2/ doF
        return self.chi2
    
   
#i do not like passing the whole std_tot
class Cluster():
    def __init__(self, event_num, std_tot, rpc, hits):
        self.event_num = event_num
        self.rpc = rpc
        self.channel_position = [-1, -1]
        self.time = [0,0]
        self.coords =  None #[x, y, z]
        self.var = None #[var_x, var_y, var_z = 0]
        self.hits = self.add_hits(hits[0], hits[1])
        self.size = max(len(self.hits[0]), len(self.hits[1]))
        self.time[1] = std_tot[self.rpc][self.channel_position[0]][self.channel_position[1]]**2
        self.coords, self.var = self.calculate_coords()
        
    def add_hits(self, phi_hits, eta_hits):
        if phi_hits:
            first_phi = min(phi_hits, key = lambda hit: hit.time)
            self.channel_position[0] = first_phi.channel
        if eta_hits:
            first_eta = min(eta_hits, key = lambda hit: hit.time)
            self.channel_position[1] = first_eta.channel
            self.time = [first_eta.time - first_eta.channel*tot_speed_factor, 0]
        return [phi_hits, eta_hits]
        
    def calculate_coords(self):
        if not self.hits[0] or not self.hits[1]:
            return None
        
        distance_per_phi_channel = 2.7625 #cm
        distance_per_eta_channel = 2.9844 #cm
        phi_channel, eta_channel = self.channel_position
        x = (phi_channel + 0.5) * distance_per_phi_channel
        y = ((31 - eta_channel) + 0.5) * distance_per_eta_channel
        var_x = (1 * distance_per_phi_channel) ** 2 / 12
        var_y = (1 * distance_per_eta_channel) ** 2 / 12
        z = RPC_heights[self.rpc]
        t = self.time
        var_t = 1 # let me ponder how to get that value here
        
        return [x, y, z], [var_x, var_y, 0]


class Reconstructor():
    def __init__(self, event_chunk, processsed_event, tolerance = None, coincidence_window = 15, tof_correction = True):
        self.event_chunk = event_chunk
        self.event_counter = 0
        self.tracks = []
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]
        self.mean_tot = [[[13 for eta in range(32)] for phi in range(64)] for rpc in range(6)]
        self.std_tot = [[[1 for eta in range(32)] for phi in range(64)] for rpc in range(6)]
        self.tof_correction = tof_correction
        self.processedEvents = processsed_event
        self.tol = [i/10 for i in range(100)] if tolerance is None else tolerance
        self.dT = []
        self.cluster_size = [[[],[]] for rpc in range(6)]
        self.collector = [[] for rpc in range(6)]
        self.chi2 = []
        self.recon = []        
        
        
        
        self.possible_reconstructions = [0 for rpc in range(6)]
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
        
    def populate_hits(self, event):
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]
        self.event_counter += 1
        skip_event = False
        for tdc in range(5):
            if event.tdcEvents[tdc].qual != 0:
                skip_event = True

            #if skip_event:
            #    continue 

            for word in event.tdcEvents[tdc].words:
                rpc, thisHit = ATools.tdcChanToRPCHit(word,tdc, self.processedEvents + self.event_counter)
                if thisHit.channel == [0]:
                    continue
                if not time_range[0] < thisHit.time < time_range[1]:
                    continue

                if thisHit.eta:
                    self.etaHits[rpc].append(thisHit)

                else:
                    self.phiHits[rpc].append(thisHit)
                                                 
    def update_event(self, event_chunk, processed_event):
        self.event_chunk = event_chunk
        self.event_counter = 0
        self.processedEvents = processed_event
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]

    def fake_speed_histogram(self, rpc_to_compare):
        fake_speed = [[] for pair in range(len(rpc_to_compare))]
        for tracks in self.tracks:
            if not tracks:
                continue
            for track in tracks:
                rpcs_detected = [rpc for rpc in range(6) if track.size[rpc] != 0]
                rpc_pos = [[] for rpc in range(6)]
                rpc_times = [0 for rpc in range(6)]
                for point in track.coordinates:
                    rpc = RPC_heights.index(point[3])
                    rpc_times[rpc] = point[0]
                    rpc_pos[rpc] = point[1:3]
               
                for idx, [rpc1, rpc2] in enumerate(rpc_to_compare):
                    if rpc1 in rpcs_detected and rpc2 in rpcs_detected:
                        fake_speed[idx].append((rpc_times[rpc2] - rpc_times[rpc1])/np.linalg.norm(np.array(rpc_pos[rpc2])-np.array(rpc_pos[rpc1])))
        return fake_speed

    def tof_histogram(self, rpc_to_compare):
        tof = [[] for pair in range(len(rpc_to_compare))]
        for tracks in self.tracks:
            if not tracks:
                continue
            for track in tracks:
                rpc_times = [0 for rpc in range(6)]
                for point in track.coordinates:
                    rpc = RPC_heights.index(point[3])
                    rpc_times[rpc] = point[0]
                rpcs_detected = [rpc for rpc in range(6) if rpc_times[rpc] != 0]
                for idx, [rpc1, rpc2] in enumerate(rpc_to_compare): #use actual for now
                    if rpc1 in rpcs_detected and rpc2 in rpcs_detected:
                        tof[idx].append(rpc_times[rpc2] - rpc_times[rpc1])
        return tof

            
    def cluster(self, time_window = 5):
        speed = 0.15204446322001586 #(25/32) ns / channel
        result = [] #[event_num, eta_time, [phi_hits, eta_hits]]
        for event_num, event in enumerate(self.event_chunk):
            self.populate_hits(event)
            event_clusters = []
            
            for rpc in range(6):
                rpc_clusters = []
                coincident_hits = [[], []]

                
                for idx, hits in enumerate([self.phiHits[rpc], self.etaHits[rpc]]):    
                    if hits:
                        temp_hits = [hits[0]]
                    else:
                        continue

                    for i in range(len(hits) - 1):    
                        if abs(hits[i+1].time - hits[i].time - speed*(hits[i+1].channel -hits[i].channel)) <= time_window: 
                            temp_hits.append(hits[i+1])
                        elif temp_hits:
                            first = min(temp_hits, key=lambda x: x.time)
                            coincident_hits[idx].append([
                                    first.time,
                                    first.channel,
                                    temp_hits
                                ])
                            temp_hits = [hits[i+1]]
                
                    if temp_hits:
                        #print([str(x) for x in temp_hits])
                        first = min(temp_hits, key=lambda x: x.time)
                        coincident_hits[idx].append([
                                    first.time,
                                    first.channel,
                                    temp_hits
                                ])
                        
                    #position clustering
                    #print("coincident_hits before", [[str(hit) for hit in temp_cluster[2]] for temp_cluster in coincident_hits[idx]])
                    for time_clusters in coincident_hits[idx]:
                        time_clusters[2] = sorted(time_clusters[2], key=lambda x: x.channel)
                        for i, hit in enumerate(time_clusters[2]):
                            if i != 0 and abs(hit.channel - previous_channel) > 1:
                                different_hits = time_clusters[2][i:].copy()
                                time_clusters[2] = time_clusters[2][:i]
                                first = min(different_hits, key=lambda x: x.time)
                                new_time_cluster = [first.time, first.channel, different_hits]
                                coincident_hits[idx].append(new_time_cluster)
                                break
                                #recursion kind of
                            previous_channel = hit.channel
                    #print("coincident_hits after", [[str(hit) for hit in temp_cluster[2]] for temp_cluster in coincident_hits[idx]])

                        
                phi_coincident, eta_coincident = coincident_hits
                leaderboard = []
                for (phi_time, phi_channel, phi_hits), (eta_time, eta_channel, eta_hits) in product(phi_coincident, eta_coincident):
                    time_difference = abs(eta_time - phi_time - self.mean_tot[rpc][phi_channel][eta_channel])/self.std_tot[rpc][phi_channel][eta_channel]
                    #is_time_coincident = abs(time_difference) < 5
                    is_time_coincident = True
                    if is_time_coincident:
                        leaderboard.append([time_difference, [phi_hits, eta_hits]])

                leaderboard = sorted(leaderboard, key=lambda x: x[0])
                original_leaderboard = leaderboard.copy()
                while leaderboard:
                    best = leaderboard.pop(0)
                    phi_hits, eta_hits = best[1]
                    rpc_clusters.append(Cluster(event_num, self.std_tot, rpc, best[1]))
                    leaderboard = [combo for combo in leaderboard if not any([hits in combo[1] for hits in [phi_hits, eta_hits]])]
                event_clusters.append(rpc_clusters)
            result.append(event_clusters)
        return result


    def reconstruct_track(self, clusters, rpc_excluded = None):
        for evt_num, event_clusters in enumerate(clusters):
            tracks = find_best_track(event_clusters, rpc_excluded)
            self.tracks.append(tracks) #inlcude event number
        return self.tracks


    #I hate this
    def reconstruct_and_extrapolate(self, dataset, chi2_region = [0, 100], only_second = False):
        # Ensure RPC is a list, even if it's a single integer
        for i, data in enumerate(dataset):
            for rpc in range(6):
                if count_entries(data) < 100:
                    E_recon = reconstruct_timed_Chi2_ByRPC(data, 3, rpc) #why is it excluded?
                    if E_recon:
                        E_recon = [E_recon] #for now take only the first
                        for track in E_recon:
                            if len(track[2]) >= 5:
                                if track[4] > chi2_region[0] and track[4] < chi2_region[1]:
                                    self.chi2.append(track[4])
                                    # self.event_of_interest.append(track)
                                    # Adding this check to see if other 5 RPCs are in reconstructed event.
                                    # This is necessary to ensure the reconstructed path is accurate.

                                    muon_coords = does_muon_hit_RPC(track[0], track[1], rpc)
                                    if muon_coords:
                                        self.possible_reconstructions[rpc] += 1
                                        self.possible_reconstructions_coords[rpc][int(muon_coords[0] / 2.7625)][int(muon_coords[1] / 2.9844)] += 1
                                        for idx, t in enumerate(self.tol):
                                            check = does_RPC_detect_muon(muon_coords, track[7], t)
                                            if check:
                                                self.successful_reconstructions[rpc][idx] += 1
                                                self.successful_reconstructed_coords[rpc][int(muon_coords[0] / 2.7625)][int(muon_coords[1] / 2.9844)] += 1
                                            else:
                                                self.failed_reconstructed_coords[rpc][int(muon_coords[0] / 2.7625)][int(muon_coords[1] / 2.9844)] += 1


    def reconstruct_and_findtof(self, dataset, rpc_comparisons):
        for i, data in enumerate(dataset):
            if count_entries(data) < 100:
                E_recon = reconstruct_timed_Chi2_ByRPC(data, 3, -1, rpc_indicies=rpc_comparisons)
                if E_recon:
                    E_recon = E_recon[0]
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
            if count_entries(filtered_event) < 100:
                result = reconstruct_timed_Chi2_modified(filtered_event, 3, max_length=max_length, exact_length=exact_length)
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
            etaHist = (scDiffs, np.array(phichannels), np.array(etachannels))
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
"""
def find_tof_time(eta, phi, slope = 0.05426554612593516, offSet = 15.8797407836404):
    if (len(set([eta.eta, phi.eta])) == 1):
        return 0
    else:
        return slope*(phi.channel-eta.channel)-offSet
"""




def clusterAFDFSD(coincident_hits):
    coincident_hits_clustered = []
    #[event_num, eta_time, [phi_hits, eta_hits]]
    for coincidence_event in coincident_hits:

        coincident_event_clustered = [coincidence_event[0], coincidence_event[1], []] #event, time, []
        hit_locations = coincidence_event[2] #hit.rpc, hit.channel, hit.time, hit.eta
        #print(hit_locations)
        phi_locations = [x for x in hit_locations if x[3] == False]
        eta_locations = [x for x in hit_locations if x[3] == True]

        phi_locations = sorted(phi_locations, key=lambda x: x[1]) #sorted by channel
        eta_locations = sorted(eta_locations, key=lambda x: x[1])

        for RPC in range(6):
            rpc_phi_clusters = []
            rpc_eta_clusters = []

            i = 0
            for index, hit in enumerate([x for x in phi_locations if x[0] == RPC]):
                if index == 0:
                    previous_element = hit[1]
                    rpc_phi_clusters.append([hit])
                else:
                    if abs(hit[1] - previous_element) > 1:
                        rpc_phi_clusters.append([hit])
                        i += 1
                    else:
                        rpc_phi_clusters[i].append(hit)
                    previous_element = hit[1]

            j = 0
            for index, hit in enumerate([x for x in eta_locations if x[0] == RPC]):
                if index == 0:
                    previous_element = hit[1]
                    rpc_eta_clusters.append([hit])
                else:
                    if abs(hit[1] - previous_element) > 1:
                        rpc_eta_clusters.append([hit])
                        j += 1
                    else:
                        rpc_eta_clusters[j].append(hit)
                    previous_element = hit[1]
            rpc_combined = [rpc_phi_clusters, rpc_eta_clusters] 
            # can they be in a different times mb??

            coincident_event_clustered[2].append(rpc_combined)

        coincident_hits_clustered.append(coincident_event_clustered)

    return coincident_hits_clustered



def extract_coords_timed_Chi2(event_clusters,max_cluster_size):

    #This function converts spatially clusters in RPCs into x and y coordinates (z given by RPC number)
    
    #Extract x and y coords of cluster in event
    distance_per_phi_channel = 2.7625 #cm
    distance_per_eta_channel = 2.9844 #cm
    
    coords = []
    for rpc, rpc_clusters in event_clusters:
        

        x_clusters = [x for x in event[2][RPC][0] if len(x)<=max_cluster_size] #phi direction
        y_clusters = [y for y in event[2][RPC][1] if len(y)<=max_cluster_size] #eta direction

        #Finding size of largest cluster, consider coordinates bad if largest cluster is larger than 6.
        # x_clusters_lengths = [len(x) for x in event[2][RPC][0]]
        # y_clusters_lengths = [len(y) for y in event[2][RPC][1]]

        # max_length = max(max(x_clusters_lengths, default=0), max(y_clusters_lengths, default=0))
        
        # if max_length > 3:
        #     return None

        x_coords = []
        y_coords = []
        
        for x_cluster in x_clusters:
            # Extract phi channels and times from the cluster
            phi_channels = [x[1] for x in x_cluster]
            phi_times = [t[2] for t in x_cluster]

            # Convert the channel number into a measurement along the RPC
            x_values = [(phi_channel + 0.5) * distance_per_phi_channel for phi_channel in phi_channels]

            # Variance in x coord
            x_var = (1 * distance_per_phi_channel) ** 2 / 12

            # Find the index of the minimum time
            min_time_index = phi_times.index(min(phi_times))

            # Append the x value corresponding to the minimum time
            x_coords.append([x_values[min_time_index], x_var, min(phi_times)])

        for y_cluster in y_clusters:
            #y_cluster = [[RPC,CHANNEL,TIME,'eta'],...]
            eta_channels_corrected = [31-y[1] for y in y_cluster] #corrected for labelling from 0 to 31.
            eta_times = [t[2] for t in y_cluster]
            y_values = [(channel_num+0.5)*distance_per_eta_channel for channel_num in eta_channels_corrected]
            
            y_var = (1*distance_per_eta_channel)**2 /12
            
            # Find the index of the minimum time
            min_time_index = eta_times.index(min(eta_times))
            
            y_coords.append([y_values[min_time_index],y_var,min(eta_times)])

        if x_coords and y_coords:

            coords.append([x_coords, y_coords])

        else: #I should think about this
            coords.append([[],[],"N"])

    #[x_coords] = [[x,err_x,x_time],...]
    
    #RPC_coords = [x_coords,y_coords]

    #coords = [[RPC1_coords],[RPC2_coords],[RPC3_coords],...]
    return coords

def generate_hit_coords_combo_Chi2(coords, RPC_heights, max_length=None, exact_length=False, combinations=None, hit_coords=None, depth=0):
    if combinations is None:
        combinations = []
    if hit_coords is None:
        hit_coords = []
    if max_length is None:
        max_length = len(coords)

    # Base case: If we've reached the end of the coords or the length condition is met
    if depth == len(coords) or len(hit_coords) == max_length:
        if not exact_length or len(hit_coords) == max_length:
            combinations.append(hit_coords.copy())
        return combinations

    # Extract x and y values for the current depth
    x_values = coords[depth][0]
    y_values = coords[depth][1]

    # If there are no x or y values at this depth, move to the next depth
    if not x_values or not y_values:
        return generate_hit_coords_combo_Chi2(coords, RPC_heights, max_length, exact_length, combinations, hit_coords, depth + 1)

    # Iterate over all combinations of x and y values
    for x in x_values:
        for y in y_values:
            if x is not None and y is not None and isinstance(x[0], (int, float)) and isinstance(y[0], (int, float)):
                hit_coords.append([x, y, RPC_heights[depth]])
                generate_hit_coords_combo_Chi2(coords, RPC_heights, max_length, exact_length, combinations, hit_coords, depth + 1)
                hit_coords.pop()

    return combinations




def extract_DT_DZ_Chi2(coords, rpc_indices):
    """
    Extracts the time difference (dT) and height difference (dZ) between specified pairs of RPCs,
    correcting for signal propagation speed along the x0 axis.

    Parameters:
    coords (list): A list of coordinates with time information for each RPC.
                   Format: [[[x0, var, time], [y0, var], z0], [[x1, var, time], [y1, var], z1], ...]
    rpc_indices (list of lists): A list of pairs of indices specifying the RPCs to compare (e.g., [[0, 1], [0, 2]]).

    Returns:
    tuple: Two lists containing dT and dZ for each pair of RPC indices.
    """
    if rpc_indices:
        dT_all = []
        dZ_all = []
        
        # Signal propagation speed correction factor
        speed_factor = 0.15204446322001586
        distance_per_phi_channel = 2.7625 #cm
        distance_per_eta_channel = 2.9844 #cm


        # Correct times for each RPC based on signal propagation speed
        times = []
        for i, coord in enumerate(coords):
            if isinstance(coord[0][2], (float, int)):
                # print(coord[1][0])
                corrected_time = coord[0][2] - ((coord[1][0] / distance_per_eta_channel) * speed_factor)
                # corrected_time = 0
                # corrected_time = coord[0][2]
                times.append((i, corrected_time))
        
        # Define heights of middle point of each RPC in cm.
        RPC_heights = [0.6, 1.8, 3.0, 61.8, 121.8, 123]

        for indices in rpc_indices:
            if len(times) > max(indices):
                t_first_rpc = times[indices[0]][1]
                t_last_rpc = times[indices[1]][1]

                dT = t_last_rpc - t_first_rpc
                dZ = RPC_heights[indices[1]] - RPC_heights[indices[0]]

                dT_all.append(dT)
                dZ_all.append(dZ)
            else:
                dT_all.append(None)
                dZ_all.append(None)
        
        return dT_all, dZ_all
    else:
        times = [[RPC,x[0][2]] for RPC, x in enumerate(coords) if isinstance(x[2], (float, int))]

        #Should already be sorted, but just in case.
        #Sort times by RPC, with RPC at lowest height at first entry.
        if len(times) > 1:
            times_sorted = sorted(times, key=lambda x: x[0])

            #print(times_sorted)

            dT = times_sorted[-1][1]-times_sorted[0][1]
            #if dT>0 this implies the particles hit the higher RPC after the lower one, so the particle is travelling upwards here.
            #Vice-versa for dT < 0 

            RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.

            first_RPC = times_sorted[0][0]
            last_RPC = times_sorted[-1][0]

            dZ = RPC_heights[last_RPC] - RPC_heights[first_RPC]
            return [dT], [dZ]
        else:
            pass


def fit_event_chi2(coordinates_with_error, rpc_indicies = None):
    #Coordinates = [[[x0,var,time],[y0,var],z0],[[x1,var,time],[y1,var],z1],...,[[x5,var,time],[y5,var],z5]]
    #Z coordinate given by height of relevant RPC.
    #Using SVD

    # Calculate dT for event, in ns
    dT, dZ = extract_DT_DZ_Chi2(coordinates_with_error, rpc_indices = rpc_indicies)
    
    coordinates = []

    for coords in coordinates_with_error:
        coordinates.append([coords[0][0],coords[1][0],coords[2]])

    centroid = np.mean(coordinates, axis=0)
    subtracted = coordinates-centroid
    # performing SVD
    try:
        _, _, V = np.linalg.svd(subtracted)
    except np.linalg.LinAlgError:
        return None, None, np.inf, None, None, None
    
    # Find the direction vector (right singular vector corresponding to the largest singular value)
    direction = V[0, :]
    # A line is defined by the average and its direction
    p0 = centroid
    d = direction

    #Work out Chi2. Minimise this to find best fit (from possible combos)

    Chi2 = 0

    i = 0 #degrees of freeedom

    for point in coordinates_with_error:
        
        i+=2
        
        z = point[2]
        x = point[0][0]
        y = point[1][0]
        x_var = point[0][1]
        y_var = point[1][1]

        z_0 = centroid[2]

        # t = (z-z_0)/d_z

        t = (z-z_0)/d[2]

        # Find expected (x,y) coordinates at that height.

        x_traj = centroid[0] + t*d[0]
        y_traj = centroid[1] + t*d[1]

        Chi2_x = (x-x_traj)**2 / x_var
        Chi2_y = (y-y_traj)**2 / y_var

        Chi2+= Chi2_x
        Chi2+= Chi2_y

    # i is number of fitted points. There are 4 fitted paramters, 2 for each x and y. 
    doF = i - 4

    Chi2 = Chi2/ doF

    return p0, d, Chi2, coordinates, dT, dZ



def fit_combination(combination):
    coordinates = []
    uncertainties = [] #no clue how to spell uncertainties
    for cluster in combination:
        #think about it
        if cluster is None:
            continue
        coordinates.append([cluster.time[0], *cluster.coords])
        uncertainties.append([cluster.time[1], *cluster.var])
    coordinates = np.array(coordinates)
    uncertainties = np.array(uncertainties)
    candidate = Track(coordinates, uncertainties)
    speed = (32/25*10**7)*np.linalg.norm(candidate.coordinates[-1][1:]-candidate.coordinates[0][1:])/(candidate.coordinates[-1][0] - candidate.coordinates[0][0])

    try:
        candidate.fit()
    except:
        return None
    return candidate

def find_best_track(event_clusters, RPC_excluded = None):
    #check empty RPCs
    empty_RPC_count = sum(1 for rpc in event_clusters if rpc == [])
    if empty_RPC_count > 3:
        return []
    #check cross chamberness
    cross_chamberness = 0
    if any(event_clusters[0:3]): #triplet
        cross_chamberness += 1
    if any(event_clusters[3:4]): #singlet
        cross_chamberness += 1
    if any(event_clusters[4:6]): #doublet
        cross_chamberness += 1

    if cross_chamberness < 2:
        return []
    #If you are testing, exclude the RPC under test
    if RPC_excluded:
        test_coords = event_clusters[RPC_excluded]
        event_clusters[RPC_excluded] = []
    #Reject big clusters
    max_cluster_size = 3
    event_clusters = [[x for x in rpc if x.size <= max_cluster_size] for rpc in event_clusters]
    for rpc in range(6):
        if event_clusters[rpc] == []:
            event_clusters[rpc] = [None]
    
    combinations = list(product(*event_clusters))
    
    #print([len(x) if x[0] else 0 for x in event_clusters])
    if all([len(x) if len(x) > 2 else 0 for x in event_clusters]):
        print([len(x) for x in event_clusters])
        print("MAYBE")

    #main loop
    possible_tracks = []
    for ind, combo in enumerate(combinations):
        if sum(1 for rpc in event_clusters if rpc == [None]) > 3:
            continue #short combination
        potential_track = fit_combination(combo)
        if potential_track and potential_track.chi2 < 10:
            possible_tracks.append(potential_track)
        #if potential_track and potential_track.chi2 < 10:
        #    possible_tracks.append(potential_track)
    #now there is a deal - which tracks am I going to return?
    #smallest chi2
    possible_tracks = sorted(possible_tracks, key = lambda x: x.chi2)
    if possible_tracks:
        return [possible_tracks[0]]
    else:
        return []
    winners = []
    while possible_tracks:
        winner = min(possible_tracks, key = lambda x: x.chi2)
        winners.append(winner)
        possible_tracks.pop(possible_tracks.index(winner))
        new_possible_tracks = []
        for ops in possible_tracks:
            if not any(taken in ops.coordinates for taken in winner.coordinates):
                new_possible_tracks.append(ops)
        if new_possible_tracks:
            print("Before:",[x.coordinates for x in possible_tracks])
            print(winner.coordinates)
            print("After:", [x.coordinates for x in new_possible_tracks])
        possible_tracks = new_possible_tracks
    return winners
    """
        if dT[-1] is not None:
            if dT[-1] > 0:
                if optimised_d[2] < 0:
                    optimised_d = np.multiply(optimised_d,-1)
            else:
                if optimised_d[2] > 0:
                    optimised_d = np.multiply(optimised_d,-1)
        if Chi2_current<max_Chi2:
            tracks.append([optimised_centroid, optimised_d, optimised_coords, combinations, Chi2_current, dT, dZ])
        else:
            return tracks
    """

def reconstruct_timed_Chi2_ByRPC(event_clusters,max_cluster_size, RPC_excluded, rpc_indicies = None):

    max_Chi2 = 6
    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]
    
    #Extract x and y coords of cluster in event
    coords = extract_coords_timed_Chi2(event, max_cluster_size)
    #[x_coords] = [[x,err_x,x_time],...]
    
    # Filter out coords of RPC under test 
    if RPC_excluded != -1:
        test_coords = event_clusters[RPC_excluded]
        event_clusters[RPC_excluded] = []

    empty_RPC_count = sum(1 for rpc in coords if rpc == [])
    if empty_RPC_count > 4:
        # print("Failed to reconstruct, not enough coords")
        return []  # Exit the function
    
    cross_chamberness = 0

    if any(event_clusters[0:3]): #triplet
        cross_chamberness += 1
    if any(event_clusters[3:4]): #singlet
        cross_chamberness += 1
    if any(event_clusters[4:6]): #doublet
        cross_chamberness += 1

    if cross_chamberness < 2:
        return []

    something = True
    tracks = []
    possible_tracks = []

    while something:
        combinations = generate_hit_coords_combo_Chi2(coords,RPC_heights)
        Chi2_current = np.inf
        optimised_coords = None
        optimised_d= None
        optimised_centroid= None
        dT = [None]

        for ind,combo in enumerate(combinations):
            print(combo)
            if len(combo) < 5:
                continue
            centroid, d, Chi2, coordinates, delta_T, delta_Z= fit_event_chi2(combo, rpc_indicies = rpc_indicies)
            if Chi2 < max_Chi2:
                possible_tracks.append([centroid, d, coordinates, combo, Chi2, delta_T, delta_Z])
            if Chi2 < Chi2_current:

                # If new fit is better than old then replace old fit properties.
                dZ = delta_Z 
                dT = delta_T
                Chi2_current = Chi2
                optimised_centroid = centroid
                optimised_d = d
                optimised_coords = coordinates
                winning_combo = combo

        #if dT>0 this implies the particles hit the higher RPC after the lower one, so the particle is travelling upwards here.
        #Vice-versa for dT < 0. The condition below make the right direction
        if dT[-1] is not None:
            if dT[-1] > 0:
                if optimised_d[2] < 0:
                    optimised_d = np.multiply(optimised_d,-1)
            else:
                if optimised_d[2] > 0:
                    optimised_d = np.multiply(optimised_d,-1)


            if Chi2_current<max_Chi2:
                # print('success')
                tracks.append([optimised_centroid, optimised_d, optimised_coords, combinations, Chi2_current, dT, dZ, test_coords])
            else:
                return tracks
        else:
            # print('Really, reallly poor Chi2 fit...')
            return tracks
        for point in winning_combo:
            try:
                rpc = RPC_heights.index(point[2])
                coords[rpc][0].remove(point[0])
                coords[rpc][1].remove(point[1])
            except:
                print("ups")
                print("Winner:", winning_combo)
                print(coords)
                print(point)    
    

def does_muon_hit_RPC(optimised_centroid, optimised_d, RPC):

    RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] 
    #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.

    # x_bar = x_centroid + d_vector * t
    # Find value of paramter t when the muon trajectory passes through the RPC height.
    
    z_0 = optimised_centroid[2]
    z = RPC_heights[RPC]

    # t = (z-z_0)/d_z

    t = (z-z_0)/optimised_d[2]

    # Find expected (x,y) coordinates at that height.

    x = optimised_centroid[0] + t*optimised_d[0]
    y = optimised_centroid[1] + t*optimised_d[1]

    # Check if these (x,y) coordinates lie within the RPC. 

    #Extract x and y coords of cluster in event
    distance_per_phi_channel = 2.7625 #cm
    distance_per_eta_channel = 2.9844 #cm

    # Max y (eta side) is 31.5 * distance_per_eta_channel
    # Max x (phi side) is 63.5 * distance_per_phi_channel

    if 0 < x < 63.5*distance_per_phi_channel and 0 < y < 31.5*distance_per_eta_channel:
        #Return coordinates where you expect the muon to hit this RPC from the reconstructed event.
        return [x,y]
    else:
        #print("Muon does not hit RPC")
        return None
    
    
def does_RPC_detect_muon(muon_coords,test_coords,tol):
    #Tolerance in units of cm. 

    #Could experiment with tolerance.
    if test_coords != [[],[],"N"]: 

        x_coords = test_coords[0]
        y_coords = test_coords[1]

        for x_set in x_coords:
            for y_set in y_coords:

                x = x_set[0]
                y = y_set[0]
    
                #If statement ensures only calculate the coords if the test_coords actually exist.

                #Offset is 2D vector that represents difference 
                offset = np.subtract(np.array([x,y]),muon_coords)

                separation = np.linalg.norm(offset)

                #print(separation)

                if separation <= tol:
                    #Say the RPC only successfully reconstructs an event 
                    #if the distance between expected hit and reconstructed hit is less than tolerance.

                    #print("RPC successfully detects hit!")
                    return separation
        
        #print("No RPC coordinates constructed pass near the expected point!")
        return False

    else:
        #print("No coordinates reconstructed by RPC")
        return False 



def reconstruct_timed_Chi2_modified(event,max_cluster_size, max_length=None, exact_length=False):

    #timed tag indicates that timing information from RPC is used to determine direction of vertical transversal of "particle" in the event.

    max_Chi2 = 10

    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]
    RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.

    #Extract x and y coords of cluster in event

    coords = extract_coords_timed_Chi2(event,max_cluster_size)

    # Count the number of empty RPCs
    empty_RPC_count = sum(1 for item in coords if item == [[], [],'N'])

    # If less than 3 elements of coords are occupied, exit the function
    if empty_RPC_count > 3:
        #print("Failed to reconstruct, not enough coords")
        return None  # Exit the function
    
    #NEED TO CHECK IF STILL CROSS CHAMBER! 

    cross_chamberness = 0

    if coords[0] != [[], [], 'N'] or coords[1] != [[], [], 'N'] or coords[2] != [[], [], 'N']:
        cross_chamberness += 1

    if coords[3] != [[], [], 'N']:
        cross_chamberness += 1

    if coords[4] != [[], [], 'N'] or coords[5] != [[], [], 'N']:
        cross_chamberness += 1

    #ITERATING OVER EVERY POSSIBLE COMBINATION OF x,y,z over all 3 RPCs (limited to one x,y per RPC).
    #Doesn't look particularly nice, but there are not many coordinates to loop over usually....

    combinations = generate_hit_coords_combo_Chi2(coords,RPC_heights, max_length=max_length, exact_length=exact_length)

    #Now for each combo in combinations, attempt to reconstruct a path. See which one gives the best trajectory.

    #If success, print parameters of fitting function.
    #If fail, print reconstruction failed.

    Chi2_current = np.inf
    optimised_coords = None
    optimised_d= None
    optimised_centroid= None
    dT = np.inf

    for ind,combo in enumerate(combinations):

        centroid, d, Chi2, coordinates, delta_T, delta_Z= fit_event_chi2(combo)
        if Chi2 < Chi2_current:

            # If new fit is better than old then replace old fit properties.
            dZ = delta_Z 
            dT = delta_T
            Chi2_current = Chi2
            optimised_centroid = centroid
            optimised_d = d
            optimised_coords = coordinates

    #if dT>0 this implies the particles hit the higher RPC after the lower one, so the particle is travelling upwards here.
    #Vice-versa for dT < 0.

    #dT = 0 case?
    
    print(dT[-1])
    if dT[-1] is not None:
        if dT[-1] > 0:
            if optimised_d[2] < 0:
                optimised_d = np.multiply(optimised_d,1)
        else:
            if optimised_d[2] > 0:
                optimised_d = np.multiply(optimised_d,0)
        print(Chi2_current)
        if Chi2_current<max_Chi2:
            return optimised_centroid, optimised_d, optimised_coords, combinations, Chi2_current, dT, dZ

        else:
            print(f"Failed to reconstruct, Chi2 too large {Chi2_current}")
            #return optimised_centroid, optimised_d, optimised_coords, combinations, residuals_current
            return None
        
        
def check_event_attributes(event,min_chamber_number,min_RPC_number):

    RPC_counter = 0
    chamber_counter = 0
    condition_1 = False
    condition_2 = False
    condition_3 = False

    for RPC in range(6):
        if RPC<3:
            if event[2][RPC][0] and event[2][RPC][1]:
                RPC_counter+=1 
                if not condition_1:
                    chamber_counter+=1
                    condition_1 = True
        elif RPC == 3:
            if event[2][RPC][0] and event[2][RPC][1]:
                RPC_counter+=1
                if not condition_2:
                    chamber_counter+=1
                    condition_2 = True
        else:
            if event[2][RPC][0] and event[2][RPC][1]:
                RPC_counter+=1
                if not condition_3:
                    chamber_counter+=1
                    condition_3 = True

    return RPC_counter >= min_RPC_number and chamber_counter >= min_chamber_number


def filter_events(events,min_chamber_number,min_RPC_number):
    filtered_events = []

    for event in events:
        if check_event_attributes(event,min_chamber_number,min_RPC_number):
            filtered_events.append(event)
    
    return filtered_events



def count_entries(lst): #why is it so complicated??
    if isinstance(lst, list):
        return sum(count_entries(sublist) for sublist in lst)
    else:
        return 1
    
    
def find_tdc_event_count(event_chunk):
    event_number = []
    for tdc in range(5):
        tot_length = 0
        for event in event_chunk:
            length = len(event.tdcEvents[tdc].words)
            tot_length += length
        event_number.append(tot_length)
    return event_number




def compile_and_plot_tof(dTs, rpc_indicies=[[0,1], [0,2], [0,3], [0,4], [0,5]], pdf_filename = "Data_output/compiled_tof_plots.pdf"):
    tof = [[] for _ in range(6)]
    for dT in dTs:
        for i in range(len(rpc_indicies)):
            tof[i].append(dT[i])
            
    def gauss(x, a, mu, sigma):
        return a * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

    
    with PdfPages(pdf_filename) as pdf:
        for i in range(len(rpc_indicies)):
            mid_average = np.mean(tof[i])
            print(f'Mid average value for RPC{i}-5: {mid_average}')

            # Perform Gaussian fit
            hist, bin_edges = np.histogram(tof[i], bins=300, range=(-20, 20))
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            try:
                popt, pcov = opt.curve_fit(gauss, bin_centers, hist, p0=[max(hist), 0, 10])

                plt.figure()
                plt.hist(tof[i], bins=1000, edgecolor='black', alpha=0.6, label='Data')
                x_fit = np.linspace(-20, 40, 1000)
                y_fit = gauss(x_fit, *popt)
                plt.plot(x_fit, y_fit, 'r-', label='Gaussian fit')
                plt.title(f'Systematically and tof Corrected, RPC{rpc_indicies[i]}')
                plt.xlim(-20, 40)
                plt.xlabel('tof / ns')
                plt.ylabel('events')
                plt.legend()

                fit_params_text = f'Amplitude = {popt[0]:.2f}\nMean = {popt[1]:.2f}\nStd Dev = {popt[2]:.2f}'
                plt.annotate(fit_params_text, xy=(0.05, 0.95), xycoords='axes fraction', verticalalignment='top')

                pdf.savefig()  # Save the current figure to the PDF
                plt.close()

                print(f'Gaussian fit parameters for RPC{rpc_indicies[i]}: amplitude = {popt[0]}, mean = {popt[1]}, std deviation = {popt[2]}')

            except RuntimeError:
                print(f'Gaussian fit failed for RPC{rpc_indicies[i]}')

        plt.figure()
        for i in range(5):
            hist, bin_edges = np.histogram(tof[i], bins=100, range=(-20, 20))
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

            try:
                popt, pcov = opt.curve_fit(gauss, bin_centers, hist, p0=[max(hist), 0, 10])
                x_fit = np.linspace(-20, 40, 1000)
                y_fit = gauss(x_fit, *popt)
                plt.plot(x_fit, y_fit, label=f'RPC{rpc_indicies[i]} fit')

            except RuntimeError:
                print(f'Gaussian fit failed for RPC{rpc_indicies[i]}')

        plt.title('Combined Gaussian Fits for RPC0-5 to RPC4-5 larger sample of 100,000')
        plt.xlim(-20, 40)
        plt.xlabel('tof / ns')
        plt.ylabel('events')
        plt.legend()

        pdf.savefig() 
        plt.close()

    return pdf_filename

def compile_and_plot_tof_chunk(dTs, rpc_indicies=[[0,1], [0,2], [0,3], [0,4], [0,5]], 
                               num_chunks=10, 
                               pdf_filename="Data_output/tof_chunks.pdf"):
    tof = [[] for _ in range(6)]
    for dT in dTs:
        for i in range(len(rpc_indicies)):
            tof[i].append(dT[i])
            
    def gauss(x, a, mu, sigma):
        return a * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

    def split_into_chunks(data, num_chunks):
        chunk_size = len(data) // num_chunks
        return [data[i*chunk_size:(i+1)*chunk_size] for i in range(num_chunks)]

    with PdfPages(pdf_filename) as pdf:
        for i in range(len(rpc_indicies)):
            tof_chunks = split_into_chunks(tof[i], num_chunks)

            for chunk_index, chunk in enumerate(tof_chunks):
                mid_average = np.mean(chunk)
                print(f'Mid average value for RPC{rpc_indicies[i]}, chunk{chunk_index + 1}: {mid_average}')

                hist, bin_edges = np.histogram(chunk, bins=300, range=(-20, 20))
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

                try:
                    popt, pcov = opt.curve_fit(gauss, bin_centers, hist, p0=[max(hist), 0, 10])

                    plt.figure(figsize=(12, 12))
                    plt.hist(chunk, bins=1000, edgecolor='black', alpha=0.6, label='Data')
                    x_fit = np.linspace(-20, 40, 1000)
                    y_fit = gauss(x_fit, *popt)
                    plt.plot(x_fit, y_fit, 'r-', label='Gaussian fit')
                    plt.title(f'Systematically and tof Corrected, RPC{rpc_indicies[i]}, chunk{chunk_index + 1}')
                    plt.xlim(-20, 40)
                    plt.xlabel('tof / ns')
                    plt.ylabel('events')
                    plt.legend()

                    fit_params_text = f'Amplitude = {popt[0]:.2f}\nMean = {popt[1]:.2f}\nStd Dev = {popt[2]:.2f}'
                    plt.annotate(fit_params_text, xy=(0.05, 0.95), xycoords='axes fraction', verticalalignment='top')

                    pdf.savefig()
                    plt.close()

                    print(f'Gaussian fit parameters for RPC{i}-{chunk_index + 1}: amplitude = {popt[0]}, mean = {popt[1]}, std deviation = {popt[2]}')

                except RuntimeError:
                    print(f'Gaussian fit failed for RPC{i}-{chunk_index + 1}')

    return pdf_filename


