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
from scipy.optimize import minimize 
import math
from itertools import groupby
from itertools import product
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

time_range = (150,350)
RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.
tot_speed_factor = 0.15204446322001586 #(25/32) ns / channel
tof_offsets = [[0, 0], [np.float64(0.25790911824685864), np.float64(0.5818488944295358)], [np.float64(0.06942782469835455), np.float64(0.5634127422432911)], [np.float64(2.17698945272121), np.float64(2.5044894876283283)], [np.float64(6.619902168271083), np.float64(2.91238011053787)], [np.float64(6.199631794023697), np.float64(2.617768586336547)]]



#time of flight, angles, chi2,...
class Track():
    def __init__(self, clusters):
        self.clusters = clusters
        self.coordinates, self.uncertainties = self.get_coordinates() 
         #t,x,y,z #var_t, var_x, var_y
        if len(self.coordinates):
            self.centroid = np.mean(self.coordinates, axis=0)
            self.cosmic = (self.coordinates[-1][0] - self.coordinates[0][0]) < 0
        self.direction = None 
        self.angles = None #phi, theta
        self.speed = 1 #c
        self.chi2 = 100 #arbitrary high value

    def __bool__(self):
        return bool(len(self.clusters))
    
    def get_coordinates(self):
        coordinates = []
        uncertainties = []
        for cluster in self.clusters:
            coordinates.append([cluster.time[0], *cluster.coords])
            uncertainties.append([cluster.time[1], *cluster.var])
        return self.tof_correction(np.array(coordinates), np.array(uncertainties))

    def tof_correction(self, coordinates, uncertainties):
        for point in coordinates:
            rpc = RPC_heights.index(point[3])
            point[0] -= tof_offsets[rpc][0]
            uncertainties[0] += tof_offsets[rpc][1]**2
        return coordinates, uncertainties
    
    def fit(self):
        try:
            _, _, V = np.linalg.svd(self.coordinates-self.centroid)
        except np.linalg.LinAlgError:
            return False
        
        self.direction = V[0, :]
        self.direction[0] = 0.0 #time should be the same
        self.direction = self.direction/np.linalg.norm(self.direction) #unit vector
        self.angles = self.find_angles()   

        if (self.coordinates[-1][0] - self.coordinates[0][0]) != 0:
            self.speed = np.linalg.norm(self.coordinates[-1][1:]-self.coordinates[0][1:])/(self.coordinates[-1][0] - self.coordinates[0][0])
            self.speed *= (32/25*10**7)/2.99792458e8
        
        self.calculate_chi2()
        return True
    
    def calculate_chi2(self):
        if self.direction[3] == 0:
            return self.chi2
            
        chi2 = 0
        i = 0 #degrees of freeedom
        for idx, point in enumerate(self.coordinates):    
            i+=3
            p = (point[3]-self.centroid[3])/self.direction[3] #parameter of the line
            for degree in range(3):
                ideal = self.centroid[degree] + p*self.direction[degree]
                real = point[degree]
                chi2+= (real-ideal)**2/self.uncertainties[idx][degree] 
        #check degrees of freedom
        doF = i - 6
        self.chi2 = chi2/ doF
        return self.chi2
    
    def find_angles(self):
        x,y,z = self.direction[1:]
        r = math.sqrt(x**2 + y**2 + z**2)
        # θ (theta) is the polar angle (from the z-axis)
        theta = math.acos(z / r)
        # φ (phi) is the azimuthal angle (from the x-axis in the xy-plane)
        phi = math.atan2(y, x)
        return (phi, theta)



class Vertex(Track): #vertex looks like a track in a triplet, double in a singlet + doublet
    def __init__(self, clusters):
        super().__init__(clusters)
        self.initial = Track([]) #initial track
        self.final = [Track([]), Track([])] #final trackS
        self.point = [] #point of intersection
        self.approach = None #distance between the two tracks
    
    def fit(self):
        initial_clusters = []
        final_clusters = [[] for rpc in range(6)] 
        #passed in different functions, so different structure
        for cluster in self.clusters:
            if cluster.rpc > 2:
                final_clusters[cluster.rpc].append(cluster)
            else:
                initial_clusters.append(cluster)
      
        self.final = find_best_track(final_clusters, 0)
        if len(self.final) < 2: # does it have two tracks?  
            return self.optimise()

        self.approach, self.point =  intersection(self.final[0], self.final[1])
      
        if self.approach > 100: #is it far away?
            return self.optimise()
        if self.point[2] > RPC_heights[3]: #is it above singlet?
            return self.optimise()
        
        
        self.initial = Track(initial_clusters)
        #include intersection
        self.initial.coordinates = np.concatenate((self.initial.coordinates, [[self.initial.centroid[0], *self.point]]), axis=0)
        self.initial.uncertainties = np.concatenate((self.initial.uncertainties, [[1, *self.initial.uncertainties[-1][1:]]]), axis=0)
        if not self.initial.fit():
            return self.optimise()
        self.calculate_chi2()
        self.initial.coordinates = self.initial.coordinates[:-1]
        self.initial.uncertainties = self.initial.uncertainties[:-1]
        self.optimise()
        return True
    
    def calculate_chi2(self):
        self.chi2 = self.initial.chi2 + self.final[0].chi2 + self.final[1].chi2
        # i think I should do something about doF but I am not sure how
        #anyway, sum will do for now
        return self.chi2
    
    def optimise(self):
        coordinates = self.coordinates[:,1:]
        if len(self.point) and self.initial:
            initial_parameters =  inverse_parametrized_vertex(self.point, self.initial.direction[1:], self.final[0].direction[1:], self.final[1].direction[1:])
        else:
            initial_parameters = [0, 0, 0, 0, 0, -0.5, -0.2, 0.2]
        f = lambda x: loss_function(coordinates, x)

        best_score = loss_function(coordinates, initial_parameters)
        best_parameters = initial_parameters
        for i in range(10):
            parameters = minimize(f, x0=initial_parameters, method="Nelder-Mead").x
            score = loss_function(coordinates, parameters)
            if score < best_score:
                best_score = score
                best_parameters = parameters
            initial_parameters = parameters + np.random.uniform(-1, 1, size=8)
        if best_score < min(50*self.chi2, 400): #400 is magic number ngl
            self.point, self.initial.direction, self.final[0].direction, self.final[1].direction = parametrized_vertex(*best_parameters)            
            for track in [self.initial, *self.final]:
                track.centroid = np.append(np.array([0.0]), self.point)
                track.direction = np.append(np.array([0.0]), track.direction)
            return True
        return False        

    def are_planer(self):
        #check if the two tracks are planar
        try:
            a = np.cross(self.final[0].direction[1:], self.final[1].direction[1:])
        except np.linalg.LinAlgError:
            return 90
        return np.arcsin(abs(np.dot(a,self.initial.direction[1:])))*180/np.pi #15 degrees deviation max

    




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
    def __init__(self, event_chunk, chunk_num):
        self.event_chunk = event_chunk
        self.event_counter = 0
        self.tracks = []
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]
        self.mean_tot = [[[13 for eta in range(32)] for phi in range(64)] for rpc in range(6)]
        self.std_tot = [[[1 for eta in range(32)] for phi in range(64)] for rpc in range(6)]
        self.chunk_num = chunk_num
        self.collector = [[] for rpc in range(6)]
        self.chi2 = []
    
        
        
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
                rpc, thisHit = ATools.tdcChanToRPCHit(word,tdc, self.chunk_num + self.event_counter)
                if thisHit.channel == [0]:
                    continue
                if not time_range[0] < thisHit.time < time_range[1]:
                    continue

                if thisHit.eta:
                    self.etaHits[rpc].append(thisHit)

                else:
                    self.phiHits[rpc].append(thisHit)
                                                 
    def update_event(self, event_chunk, chunk_num):
        self.event_chunk = event_chunk
        self.event_counter = 0
        self.chunk_num = chunk_num
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]

    def fake_speed_histogram(self, rpc_to_compare):
        fake_speed = [[] for pair in range(len(rpc_to_compare))]
        for tracks in self.tracks:
            if not tracks:
                continue
            for track in tracks:
                rpc_pos = [[] for rpc in range(6)]
                rpc_times = [0 for rpc in range(6)]
                for point in track.coordinates:
                    rpc = RPC_heights.index(point[3])
                    rpc_times[rpc] = point[0]
                    rpc_pos[rpc] = point[1:3]
                rpcs_detected = [rpc for rpc in range(6) if rpc_times[rpc] != 0]

                for idx, [rpc1, rpc2] in enumerate(rpc_to_compare):
                    if rpc1 in rpcs_detected and rpc2 in rpcs_detected:
                        value = (rpc_times[rpc2] - rpc_times[rpc1])/np.linalg.norm(rpc_pos[rpc2]-rpc_pos[rpc1])
                        if value is not None and np.isfinite(value):
                            fake_speed[idx].append(value)
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

                #for both eta and phi
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


    def reconstruct_tracks(self, chunk, chunk_num):
        self.update_event(chunk, chunk_num)
        clusters = self.cluster()
        for evt_num, event_clusters in enumerate(clusters):
            tracks = find_best_track(event_clusters, (self.chunk_num, evt_num))
            self.tracks.append(tracks) #inlcude event number
        return self.tracks
    #polar, azimuthal, solid angle, [-pi,pi]                    
                

def find_best_track(event_clusters, event_id, RPC_excluded = None):
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
    event_clusters = [[cluster for cluster in rpc if cluster.size <= max_cluster_size] for rpc in event_clusters]
    for rpc in range(6):
        if event_clusters[rpc] == []:
            event_clusters[rpc] = [None]
            #None must be included to generate combos
    combinations = list(product(*event_clusters))
    
    #very easy way how to detect candidates for vertex
    cluster_numbers = [len(x) if x != [None] else 0 for x in event_clusters]
    if cluster_numbers[:3] == [1,1,1] and all([cl > 1 for cl in cluster_numbers[3:]]):
        vertex = Vertex([cluster for clusters in event_clusters for cluster in clusters])
        if vertex.fit():
            return [vertex]
    #main loop
    possible_tracks = []
    for ind, combo in enumerate(combinations):
        potential_track = Track([cluster for cluster in combo if cluster != None]) #combo = [None, cluster1, ...] len = 6
        if potential_track.fit() and potential_track.chi2 < 15:
            possible_tracks.append(potential_track)

    possible_tracks = sorted(possible_tracks, key = lambda x: x.chi2)       
    #now there is a deal - which tracks am I going to return?
    #smallest chi2
    winners = []
    while possible_tracks:
        winner = possible_tracks.pop(0)
        winners.append(winner)
        new_possible_tracks = []
        for ops in possible_tracks:
            ops_coordinates = ops.coordinates.tolist()
            winner_coordinates = winner.coordinates.tolist()
            if not any(taken in ops_coordinates for taken in winner_coordinates):
                new_possible_tracks.append(ops)
        possible_tracks = new_possible_tracks
    return winners
            


def intersection(track_a, track_b): #find the closest distance and point of the intersection
    a = track_a.direction[1:]
    b = track_b.direction[1:] #direction vectors
    ca = track_a.centroid[1:]
    cb = track_b.centroid[1:]
    if np.linalg.norm(np.cross(a,b)):
        normal = np.cross(a, b)/np.linalg.norm(np.cross(a,b)) #looking from b
    else:
        return 10000.0, [0.0, 0.0, 0.0]
    closest_distance = np.dot(normal, ca - cb) #conserve sign
    a_parameter = np.cross(cb - ca + normal*closest_distance, b)/np.cross(a,b) #parameters look good
    closest_point = a*a_parameter[0]+ca-normal*closest_distance/2
    return abs(closest_distance), closest_point

def normalize(vector):
    return vector/np.linalg.norm(vector)

def sigmoid(x):
    return 1/(1+np.exp(-x))

def parametrized_vertex(p1,p2,p3, theta, phi, gamma, a, b):
    """
    vertex in the box - 3 parameters
    plane angles - 2 parameters
    you need to choose the opening angles - 3 parametes
    8 tof
    """
    distance_per_phi_channel = 2.7625 #cm
    distance_per_eta_channel = 2.9844 #cm
    dimensions = [distance_per_phi_channel*64, distance_per_eta_channel*32, 123] #cm
    
    vertex_point = np.array([sigmoid(p1)*dimensions[0], sigmoid(p2)*dimensions[1], (0.1+0.5*sigmoid(p3))*dimensions[2]]) #ready
    theta, phi = 0.4*np.arctan(theta), 2*np.arctan(phi) #
    n = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), np.sin(theta)]) #normal vector
    gamma, a, b = np.arctan(theta), np.arctan(a),np.arctan(b) #between -pi/2 < arctan(a) < pi/2 
    #find the vertex point
    ezp = normalize(np.array([0,0,1]) - np.dot(np.array([0,0,1]),n)*n)
    eyp = np.cross(n, ezp)
    #find the direction
    initial = -ezp*np.cos(gamma) - eyp*np.sin(gamma)
    finala = np.cos(gamma+a)*ezp + np.sin(gamma+a)*eyp
    finalb = np.cos(gamma+b)*ezp + np.sin(gamma+b)*eyp
    return  vertex_point, initial, finala, finalb #initial, finala, finalb

def inverse_parametrized_vertex(vertex_point, initial, finala, finalb):
    # Known constants
    distance_per_phi_channel = 2.7625  # cm
    distance_per_eta_channel = 2.9844  # cm
    dimensions = [distance_per_phi_channel * 64, distance_per_eta_channel * 32, 123]  # cm
    
    # Step 1: Recover p1, p2, p3 using the inverse of the sigmoid
    def sigmoid_inv(x):
        return np.log(x / (1 - x))
    
    p1 = sigmoid_inv(vertex_point[0]  / (dimensions[0] ))
    p2 = sigmoid_inv(vertex_point[1] / dimensions[1] )
    p3 = sigmoid_inv(vertex_point[2] / dimensions[2] )

    # Step 2: Recover theta and phi from the normal vector n
    # The normal vector n can be derived from initial or finala/finalb. For now, assuming `initial`.
    n = normalize(np.cross(initial, finala))  # Normal vector n from initial and z-axis
    
    theta = np.arcsin(n[2])  # theta from n_z = sin(theta)
    phi = np.arctan2(n[1], n[0])  # phi from n_x and n_y = cos(theta)
    
    # In the original function: theta = 0.2 * np.arctan(theta), phi = 2 * np.arctan(phi)
    theta = np.tan(theta) / 0.4  # Recover original theta
    phi = np.tan(phi) / 2  # Recover original phi
    
    # Step 3: Recover gamma, a, and b from initial, finala, finalb
    gamma = np.arctan2(-initial[1], -initial[0])  # Use initial to recover gamma
    
    # Using finala and finalb to recover a and b
    a = np.arctan2(finala[1], finala[0]) - gamma  # Recover a
    b = np.arctan2(finalb[1], finalb[0]) - gamma  # Recover b
    
    return p1, p2, p3, theta, phi, gamma, a, b


def loss_function(points, parameters):
    """
    find intersection with RPCs, find the distance to the point, chi2
    passing points without time
    """
    distance_per_phi_channel = 2.7625 #cm
    distance_per_eta_channel = 2.9844 #cm
    var_x = (1 * distance_per_phi_channel) ** 2 / 12
    var_y = (1 * distance_per_eta_channel) ** 2 / 12
        
    vertex, initial, final_a, final_b = parametrized_vertex(*parameters)
    score = 0
    for point in points:
        height = point[2]
        if height < vertex[2]: #check initial
            t = (point[2] - vertex[2])/initial[2]    
            x = vertex[0] + t*initial[0]
            y = vertex[1] + t*initial[1]
            score += (point[0] - x)**2/var_x + (point[1] - y)**2/var_y
        else: #check final
            t = (point[2] - vertex[2])/final_a[2]    
            x = vertex[0] + t*final_a[0]
            y = vertex[1] + t*final_a[1]
            score_a = (point[0] - x)**2/var_x + (point[1] - y)**2/var_y
            t = (point[2] - vertex[2])/final_b[2]    
            x = vertex[0] + t*final_b[0]
            y = vertex[1] + t*final_b[1]
            score_b = (point[0] - x)**2/var_x + (point[1] - y)**2/var_y 
            score += min(score_a, score_b)
    return score       


    