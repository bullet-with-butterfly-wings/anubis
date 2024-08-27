import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # or 'Qt5Agg', 'GTK3Agg', etc.
import mplhep as hep
hep.style.use([hep.style.ATLAS])
import sys
# import ANUBIS_triggered_functions as ANT
import matplotlib.backends.backend_pdf
import numpy as np
# from scipy.stats import normpip install pillow
sys.path.insert(1, 'Osiris Temp\processing\python')
from itertools import groupby
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def find_tof_time(eta, phi, slope = 0.05426554612593516, offSet = 15.8797407836404):
    if (len(set([eta.eta, phi.eta])) == 1):
        return 0
    else:
        return slope*(phi.channel-eta.channel)-offSet

def FindCoincidentHits(etaHits, phiHits, time_window, tof_correction = True, slope = 0.05426554612593516, offSet = 15.8797407836404):
    channels = []

    for RPC in range(6):
        channels += [hit for hit in etaHits[RPC] if 150 <= hit.time <= 300 and hit.channel != 0]
        
        channels += [hit for hit in phiHits[RPC] if 150 <= hit.time <= 300 and hit.channel != 0]
        
    event_sorted = sorted(channels, key=lambda rpcHit: (rpcHit.event_num, rpcHit.time))
    grouped_and_sorted = {key: list(group) 
                          for key, group in groupby(event_sorted, lambda rpcHit: rpcHit.event_num)}
    #print("Grouped")
    #p
    #rint(grouped_and_sorted)
    coincident_hits = []

    for event_num, hits in grouped_and_sorted.items():
        temp_hits = []

        for i in range(len(hits) - 1):
            if tof_correction:
                if hits[i].eta and not hits[i+1].eta:
                    correction = find_tof_time(hits[i], hits[i+1], slope, offSet)
                elif not hits[i] and hits[i+1].eta:
                    correction = find_tof_time(hits[i+1], hits[i], slope, offSet)
                else:
                    correction = 0
            else:
                correction = 0
            
            if abs(hits[i+1].time - hits[i].time - correction) <= time_window: #can it be out of order due to the correction?
                temp_hits.append(hits[i])
                temp_hits.append(hits[i+1])
            #I think this is wrong?? It can be both phi, 
            # Does not represent quality Better??
        #shift
        if temp_hits:
            unique_hits = { (hit.channel, hit.time, hit.eta, hit.event_num, hit.rpc): hit for hit in temp_hits }.values()
            eta_hits = [hit for hit in unique_hits if hit.eta]
            if eta_hits:
                time_bin = min(hit.time for hit in eta_hits) 
                #highly speculative
            else:
                time_bin = 1 #again
        
            coincident_hits.append([
                event_num,
                time_bin,
                [[hit.rpc, hit.channel, hit.time, hit.eta] for hit in unique_hits]
            ])
    #print("Coin")
    #print(coincident_hits)
    return coincident_hits


def cluster(coincident_hits):
    coincident_hits_clustered = []

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



def extract_coords_timed_Chi2(event,max_cluster_size):

    #This function converts spatially clusters in RPCs into x and y coordinates (z given by RPC number)
    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]

    #Extract x and y coords of cluster in event
    distance_per_phi_channel = 2.7625 #cm
    distance_per_eta_channel = 2.9844 #cm
    
    coords = []
    for RPC in range(6):
        
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

        else:
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

    i = 0 

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


def reconstruct_timed_Chi2_ByRPC(event,max_cluster_size, RPC_excluded, rpc_indicies = None):

    max_Chi2 = 10

    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]
    RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.


    #Extract x and y coords of cluster in event
    coords = extract_coords_timed_Chi2(event,max_cluster_size)

    test_coords = -1

    # Filter out coords of RPC under test 
    if RPC_excluded != -1:
        test_coords = coords[RPC_excluded]

        coords[RPC_excluded] = [[],[],"N"] 
    # print(coords)
    # Count the number of empty RPCs
    empty_RPC_count = sum(1 for item in coords if item == [[], [],'N'])

    # If less than 3 elements of coords are occupied, exit the function
    if empty_RPC_count > 4:
        # print("Failed to reconstruct, not enough coords")
        return None  # Exit the function
    
    cross_chamberness = 0

    if coords[0] != [[], [], 'N'] or coords[1] != [[], [], 'N'] or coords[2] != [[], [], 'N']:
        cross_chamberness += 1

    if coords[3] != [[], [], 'N']:
        cross_chamberness += 1

    if coords[4] != [[], [], 'N'] or coords[5] != [[], [], 'N']:
        cross_chamberness += 1

    if cross_chamberness < 2:
        # print("Failed to reconstruct, too few chambers")
        return None

    #ITERATING OVER EVERY POSSIBLE COMBINATION OF x,y,z over all 3 RPCs (limited to one x,y per RPC).
    #Doesn't look particularly nice, but there are not many coordinates to loop over usually....

    combinations = generate_hit_coords_combo_Chi2(coords,RPC_heights)

    if len(combinations) > 20:
        return None
    #Now for each combo in combinations, attempt to reconstruct a path. See which one gives the best trajectory.

    #If success, print parameters of fitting function.
    #If fail, print reconstruction failed.

    Chi2_current = np.inf
    optimised_coords = None
    optimised_d= None
    optimised_centroid= None
    dT = [None]

    for ind,combo in enumerate(combinations):
        if len(combo) < 5:
            continue
        centroid, d, Chi2, coordinates, delta_T, delta_Z= fit_event_chi2(combo, rpc_indicies = rpc_indicies)
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
    if dT[-1] is not None:
        if dT[-1] > 0:
            if optimised_d[2] < 0:
                optimised_d = np.multiply(optimised_d,-1)
        else:
            if optimised_d[2] > 0:
                optimised_d = np.multiply(optimised_d,-1)


        if Chi2_current<max_Chi2:
            # print('success')
            return optimised_centroid, optimised_d, optimised_coords, combinations, Chi2_current, dT, dZ, test_coords

    else:
        # print('Really, reallly poor Chi2 fit...')
        return None
    
    
    

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
    
    if dT[-1] is not None:
        if dT[-1] > 0:
            if optimised_d[2] < 0:
                optimised_d = np.multiply(optimised_d,1)
        else:
            if optimised_d[2] > 0:
                optimised_d = np.multiply(optimised_d,0)

        if Chi2_current<max_Chi2:
            return optimised_centroid, optimised_d, optimised_coords, combinations, Chi2_current, dT, dZ

        else:
            #print("Failed to reconstruct, Chi2 too large")
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
