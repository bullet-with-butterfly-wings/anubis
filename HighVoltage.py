#!/usr/bin/env python
# coding: utf-8

# In[11]:


import overview
import importlib
import rawFileReader
import pickle
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

dir_path = "/usera/jd2052/Documents/anubis/" # insert your directory path


# In[8]:


#Save chunks
files = ["proAnubis_240815_1759.raw", "proAnubis_240818_1125.raw","proAnubis_240818_1125.raw", "proAnubis_240818_1325.raw", "proAnubis_240829_2026.raw", "proAnubis_240829_2026.raw"]

rpc = 4
storage_name = f"chunks_hv{rpc}.pkl"
file_name = files[rpc]

hv_file = 'data/hvScan.csv'  # Replace with your file path
hv_data = pd.read_csv(hv_file, usecols=["start_"+str(rpc),"end_"+str(rpc),"voltage_"+str(rpc)])

#rpc 0: 2024-08-15 17:24:22	
#rpc 1	2024-08-18 11:39:42
#rpc 2	2024-08-18 12:22:21 


# In[10]:


fReader = rawFileReader.fileReader(dir_path+"data//"+file_name) # load in the classs object
total_chunks = []
for start, end, voltage in hv_data.values:
    if type(start) == float:
        break
    print("Voltage:", voltage)
    start = "2024-08-29 " + start 
    end = "2024-08-29 "+ end
    chunks, times, fReader = overview.get_chunks(file_name, fReader=fReader, start=start, end=end)
    total_chunks.append((voltage, times, chunks))


# In[23]:


with open(storage_name, "wb") as outp:
    pickle.dump(total_chunks, outp)
print("Chunks Saved")


# In[ ]:


with open(storage_name, "rb") as inp:
    total_chunks = pickle.load(inp)
print("Chunks Loaded")


# In[11]:



complete_data = []

with PdfPages(f"hit_time_histograms{rpc}.pdf") as pdf: #trash pdf, so I do not need to click on
    for voltage, times, chunks in total_chunks:
        results_dict = {}
        print(f"# Chunks for {round(voltage)}:", len(chunks))
        if len(chunks) == 1:
            print("Only one chunk, skipping")
            continue
        
        cluster_size, error, hist = overview.cluster_size(chunks, residual=False)
        results_dict["cluster_size"] = (cluster_size, error, hist)
        results_dict["efficiency"] = overview.efficiency(chunks, residual = False, pdf = pdf)
        
        
        good, bad = overview.abs_bvg_hits(chunks[:100], times[:100], per_rpc=True, pdf = pdf)[rpc]
        
        results_dict["counts"] = (good, bad)
        
        hist, std = overview.hit_time_hist(chunks, per_rpc= True, pdf=pdf)
        results_dict["hit_time_hist"] = (hist, std)
        results_dict["voltage"] = voltage
        if results_dict:
            complete_data.append(results_dict)


        


# In[15]:


print(complete_data[-1])


# In[12]:

"""
import matplotlib.pyplot as plt
import numpy as np
selected_voltages = [5600, 5800, 6000]
for volt_data in complete_data:
    if volt_data["voltage"] in selected_voltages and "cluster_size" in volt_data.keys():
        hist, bins = np.histogram(volt_data["cluster_size"][2][rpc], bins=range(1, max(volt_data["cluster_size"][2][rpc])+2))
        plt.step(bins[:-1], hist, label="Voltage: "+volt_data["voltage"]+"V", alpha=0.5)
plt.legend()
plt.xlabel("Cluster Size")
plt.ylabel("Counts")
plt.yscale("log")
plt.title(f"Cluster Size Distribution RPC {rpc}")
plt.show()
"""


# In[16]:


import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))  # Create a new figure

# Plot In-peak and Off-peak hits in the first figure

good = []
bad = []
voltage = []
for volt_data in complete_data:
    if "counts" not in volt_data.keys():
        continue
    good.append(sum(volt_data["counts"][0])/100)
    bad.append(sum(volt_data["counts"][1])/100)
    voltage.append(volt_data["voltage"])

plt.plot(voltage, good, label="In-peak hits", color="blue")
plt.plot(voltage, bad, label="Off-peak hits", color="red")
plt.title(f"Off-peak and In-peak hits (RPC {rpc})")
plt.xlabel("Voltage (V)")
plt.ylabel("Counts per chunk")
plt.legend()
plt.show()

print("Good",good )
print("Bad", bad)
print("Voltage", voltage)


# In[109]:


# Plot Efficiency in the third figure
"""
plt.figure(figsize=(10, 6))  # Create a new figure
efficiencies = []
voltage = []
for volt_data in complete_data:
    if "efficiency" not in volt_data.keys():
        continue
    efficiencies.append(volt_data["efficiency"])
    voltage.append(volt_data["voltage"])
plt.plot(voltage, efficiencies, label="Efficiency", color="red")
plt.title(f"Efficiency (RPC {rpc})")
plt.xlabel("Voltage (V)")
plt.ylabel("Efficiency")
plt.legend()
plt.show()

print("Efficiencies:", efficiencies)
print("Voltages:", voltage)
"""

# In[9]:


# Plot Cluster size in the second figure
"""
clusters = []
errors = []
voltage = []
for volt_data in complete_data:
    if "cluster_size" not in volt_data.keys():
        continue
    clusters.append(volt_data["cluster_size"][0][rpc])
    errors.append(volt_data["cluster_size"][1][rpc])
    voltage.append(volt_data["voltage"])

plt.figure(figsize=(10, 6))  # Create a new figure
plt.plot(voltage, clusters, label="Cluster size", color="red")
plt.errorbar(voltage, clusters, yerr=errors, fmt='o', color="red")
plt.title(f"Cluster size (RPC {rpc})")
plt.xlabel("Voltage (V)")
plt.ylabel("Cluster size")
plt.legend()
plt.show()

print("Cluster size:", clusters)
print("Errors:", errors)
print("Voltages", voltage)


# In[22]:


#Plot Histograms
rpc = 4
selected_voltages = [5600, 5800, 6000]
plt.figure(figsize=(10, 6))  # Create a new figure
for i, volt_data in enumerate(complete_data):
    if volt_data["voltage"] in selected_voltages and "hit_time_hist" in volt_data.keys():
        hist = volt_data["hit_time_hist"][0][0]
        plt.plot(hist[rpc][0][:-1], hist[rpc][1], label=f"Voltage {volt_data['voltage']} V")
        plt.title(f"Hit time histogram (RPC {rpc})")
        plt.xlabel("Time (25/32 ns)")
        plt.yscale("log")
        plt.ylabel("Counts (normalised)")
plt.legend()
plt.show()
"""

# In[40]:


stds = []
voltages = []

for volt_data in complete_data:
    if "hit_time_hist" in volt_data.keys():
        stds.append(volt_data["hit_time_hist"][1][rpc])
        voltages.append(volt_data["voltage"])

plt.plot(voltages, stds , label="Standard deviation")
plt.xlabel("Voltage (V)")
plt.ylabel("Standard deviation")
plt.title(f"Spread of the hit time peak (RPC {rpc})")
plt.legend()
plt.show()

print("Stds:", stds)
print("Voltages:", voltages)


# In[46]:

"""
#Save chunks

files = ["proAnubis_240815_1759.raw", "proAnubis_240818_1125.raw", "proAnubis_240818_1325.raw"]


storage_name = f"chunks_hv{rpc}.pkl"
file_name = files[rpc]

hv_file = 'data/hvScan.csv'  # Replace with your file path
hv_data = pd.read_csv(hv_file)
#hv_data = pd.read_csv(hv_file, usecols=["start_"+str(rpc),"end_"+str(rpc),"voltage_"+str(rpc)])

#rpc 0: 2024-08-15 17:24:22	
#rpc 1	2024-08-18 11:39:42
#rpc 2	2024-08-18 12:22:21 


import numpy as np
import matplotlib.pyplot as plt
cl0 = [np.float64(1.0339622641509434), np.float64(1.045432730763078), np.float64(1.0991290983606556), np.float64(1.1098737353933026), np.float64(1.158373989173294), np.float64(1.199706026457619), np.float64(1.243627828510654), np.float64(1.3051052260598126), np.float64(1.3432812085873311), np.float64(1.4434467018370936), np.float64(1.5500040730985418), np.float64(1.6012395661518468), np.float64(1.7499497487437186), np.float64(1.9432747311988978), np.float64(2.1125484678996673), np.float64(2.418690392760291), np.float64(2.816230601735043), np.float64(3.3619210977701544), np.float64(4.264796643667002), np.float64(5.270159270760293)]
err0 = [np.float64(0.011126854946426285), np.float64(0.00797484375042157), np.float64(0.010639097472394752), np.float64(0.006399245040924702), np.float64(0.005373077722241347), np.float64(0.011612885000507504), np.float64(0.011525592304870512), np.float64(0.012537685498550555), np.float64(0.011259031916725452), np.float64(0.01603256469828638), np.float64(0.016418504051676275), np.float64(0.016081504890198023), np.float64(0.019236149879349896), np.float64(0.02214773563960034), np.float64(0.026127083422149457), np.float64(0.028637597625313208), np.float64(0.03627862360495955), np.float64(0.03986771546091226), np.float64(0.05496790069215523), np.float64(0.06452445127451842)]
eff_0 = [0.7142857142857143, 0.631578947368421, 0.7407407407407407, 0.800498753117207, 0.837734404898584, 0.8513661202185793, 0.8619281045751634, 0.8058705803869246, 0.7712418300653595, 0.7727272727272727, 0.8007478632478633, 0.8371369294605809, 0.06060606060606061, 0.8083756345177665, 0.7622107969151671, 0.6454965357967667, 0.6234747239976758, 0.5983076157292185, 0.5689544579858884, 0.6159267089499648]
good_0 = [27.23, 176.75, 187.79, 208.36, 217.25, 219.96, 231.34, 237.86, 244.11, 260.15, 282.66, 307.92, 334.27, 358.79, 387.47, 420.82, 458.28, 558.12, 707.28, 952.21]
bad_0 = [0.36, 7.19, 10.49, 13.81, 15.59, 17.85, 20.35, 25.15, 26.74, 27.65, 39.5, 38.18, 45.49, 48.32, 61.55, 60.18, 68.13, 77.08, 101.84, 110.29]
std_0 = [np.float64(116.44640183214655), np.float64(147.14398853196343), np.float64(130.4686001108644), np.float64(60.85581179078705), np.float64(147.10844128432203), np.float64(175.56184808696304), np.float64(195.82146789505632), np.float64(203.24078063739753), np.float64(222.7283138879541), np.float64(216.74833196589313), np.float64(233.39542971654942), np.float64(238.1559128032108), np.float64(248.41570803066662), np.float64(254.446597678295), np.float64(254.71945322400353), np.float64(264.04218057542863), np.float64(254.63390691215753), np.float64(260.4522079461875), np.float64(261.20382309648596), np.float64(268.37842991368245), np.float64(265.0101410268362), np.float64(256.88603378265293), np.float64(246.9555944252669)]
voltages_0 = [3000.0, 4000.0, 4500.0, 4750.0, 5000.0, 5100.0, 5200.0, 5250.0, 5300.0, 5350.0, 5400.0, 5450.0, 5500.0, 5550.0, 5600.0, 5650.0, 5700.0, 5750.0, 5800.0, 5850.0, 5900.0, 5950.0, 6000.0]


cl1 = [np.float64(1.2212389380530972), np.float64(1.1514522821576763), np.float64(1.047808764940239), np.float64(1.074378831881252), np.float64(1.0793083854208945), np.float64(1.130030772943945), np.float64(1.1806598517513687), np.float64(1.2741918045194076), np.float64(1.3286349856293556), np.float64(1.3650258038904328), np.float64(1.4437348498974454), np.float64(1.5276988398916787), np.float64(1.6550848547184367), np.float64(1.7486767031002957), np.float64(1.9628453593196833), np.float64(2.2378623115171603), np.float64(2.5977705099267583), np.float64(3.288475477168398), np.float64(4.1703055084091964), np.float64(5.144818364163138)]
err1 = [np.float64(0.15365365868781636), np.float64(0.050306133269046746), np.float64(0.01085851940062489), np.float64(0.012610236466712046), np.float64(0.007588714747234421), np.float64(0.009364968167379703), np.float64(0.010283760210871348), np.float64(0.012014437176205687), np.float64(0.013788835005604861), np.float64(0.012602875350185448), np.float64(0.017094625432900766), np.float64(0.016288766736714727), np.float64(0.01812416969437741), np.float64(0.019781668576846993), np.float64(0.023790335485086862), np.float64(0.029631820011282582), np.float64(0.03495556082164944), np.float64(0.044105284020961893), np.float64(0.05386122345286682), np.float64(0.04904885980725799)]
eff_1 = [0, 0.6363636363636364, 0.6, 0.6492890995260664, 0.6837416481069042, 0.7642526964560863, 0.8214285714285714, 0.7619439868204283, 0.7778810408921933, 0.8096013018714402, 0.8154020385050963, 0.8282910874897792, 0.7958030669895076, 0.7861271676300579, 0.7273449920508744, 0.6556451612903226, 0.6175889328063241, 0.6141078838174274, 0.6260543580131209, 0.6921641791044776]
good_1 = [23.11, 101.19, 143.96, 161.91, 178.94, 195.78, 209.24, 210.39, 224.92, 246.56, 265.62, 278.61, 294.74, 311.62, 325.12, 350.22, 413.4, 535.85, 680.4, 946.93]
bad_1 = [0.82, 2.85, 4.97, 10.8, 11.69, 13.95, 19.5, 22.54, 28.8, 30.71, 39.49, 39.49, 43.2, 48.55, 57.97, 48.03, 68.53, 78.68, 94.05, 143.15]
std_1 = [np.float64(33.62848900345012), np.float64(75.40251771826425), np.float64(20.23976986093638), np.float64(98.9403675674452), np.float64(157.65312727971587), np.float64(25.425914261172096), np.float64(155.12154453706253), np.float64(122.58264578206251), np.float64(143.98342398843823), np.float64(176.8360062555137), np.float64(190.7766510634326), np.float64(208.47690585506578), np.float64(231.04675647537343), np.float64(239.5889561940645), np.float64(244.3349516634308), np.float64(242.81730853509387), np.float64(255.74022309532998), np.float64(259.7230028227785), np.float64(263.26304771075047), np.float64(263.73677604056456), np.float64(272.60875767022014), np.float64(264.79901299335796), np.float64(270.115786313314), np.float64(267.52328190045745), np.float64(259.5640546638369), np.float64(256.75948891329324)]
voltages_1 = [4000, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850, 5900, 5950, 6000]

cl3 = [np.float64(1.3109243697478992), np.float64(1.4304932735426008), np.float64(1.088797108931337), np.float64(1.0732620320855615), np.float64(1.0889143293264645), np.float64(1.140652965500216), np.float64(1.2034606205250598), np.float64(1.300078919270229), np.float64(1.3564630027632791), np.float64(1.4238995328104238), np.float64(1.585788485880145), np.float64(1.614696409399976), np.float64(1.7720491102649778), np.float64(1.9402696793002916), np.float64(2.200837122032186), np.float64(2.4071046264388887), np.float64(2.8758147771372595), np.float64(3.5219217537402994), np.float64(4.607378273061656), np.float64(5.997637872773203)]
err3 = [np.float64(0.13743948558122332), np.float64(0.17042343068733465), np.float64(0.029088649923025792), np.float64(0.014510821401138367), np.float64(0.01037693249897322), np.float64(0.00932813750248753), np.float64(0.01002816662956866), np.float64(0.012962796812449724), np.float64(0.015681357316640058), np.float64(0.015622954630547714), np.float64(0.021914920916389528), np.float64(0.019423017007949815), np.float64(0.02219151199022226), np.float64(0.024554762448446458), np.float64(0.028902655790087), np.float64(0.031205380199562274), np.float64(0.03841785594435168), np.float64(0.047352697168423274), np.float64(0.057293592625362544), np.float64(0.06856915514193454)]
eff_3 = [0.3333333333333333, 0, 0.6046511627906976, 0.6592592592592592, 0.7282229965156795, 0.8059701492537313, 0.8470209339774557, 0.8125894134477826, 0.8591352859135286, 0.7796610169491526, 0.7790262172284644, 0.83729662077597, 0.8648960739030023, 0.8652173913043478, 0.8177676537585421, 0.7995712754555199, 0.7423913043478261, 0.6504524886877828, 0.6239035087719298, 0.6424180327868853]
good_3 = [39.41, 118.95, 167.0, 178.56, 195.79, 211.84, 229.26, 247.03, 251.18, 267.39, 281.43, 314.74, 349.06, 390.47, 420.62, 477.23, 521.25, 645.4, 786.85, 1058.32]
bad_3 =[1.24, 3.28, 5.57, 9.0, 12.72, 14.29, 17.94, 22.83, 22.58, 30.95, 27.29, 30.73, 42.35, 41.41, 49.63, 57.62, 68.54, 78.72, 95.06, 103.75]
std_3 = [np.float64(22.580777820499762), np.float64(36.26761181903454), np.float64(115.76890863771862), np.float64(45.15559029300194), np.float64(15.544845883103225), np.float64(50.23009591777843), np.float64(109.6711771245793), np.float64(117.64960279069595), np.float64(139.9139886519241), np.float64(164.5487685457531), np.float64(179.88285706818914), np.float64(197.46726684488053), np.float64(204.03595753014383), np.float64(213.83843936997752), np.float64(222.52963942577753), np.float64(228.13897591693993), np.float64(236.74615068110663), np.float64(236.3079104812058), np.float64(241.65472667069506), np.float64(252.93947698782003), np.float64(242.48333735286218), np.float64(245.4035123716854), np.float64(241.87816539584296), np.float64(240.27277861641534), np.float64(230.55177350620454), np.float64(225.56577025801803)]
voltages_3 = [4000, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850, 5900, 5950, 6000]

good_4 = [3.28, 18.78, 36.45, 45.05, 45.68, 81.54, 120.08, 206.02, 232.52, 228.09, 234.95, 286.25, 257.07, 409.36, 378.03, 422.63, 549.94, 543.93, 541.65, 707.08]
bad_4 = [0.41, 1.22, 2.91, 3.45, 3.74, 7.42, 12.18, 23.86, 33.13, 32.67, 31.82, 46.7, 36.0, 64.93, 68.7, 72.77, 103.94, 87.74, 82.86, 107.52]
voltage_4 = [4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5450, 5500, 5550, 5600, 5650, 5700, 5750, 5800, 5850, 5900, 5950, 6000]



cluster = [cl0, cl1, cl3]
error = [err0, err1, err3]
eff = [eff_0, eff_1, eff_3]
good = [good_0, good_1, good_3]
bad = [bad_0, bad_1, bad_3]
stds = [std_0, std_1, std_3]
voltages = [voltages_0, voltages_1, voltages_3]

#plot cluster size

colour = ["red", "blue", "green"]
labels = ["bottom triplet", "middle triplet", "singlet"]
plt.figure(figsize=(10, 6))  # Create a new figure
for rpc in range(3):
    print(rpc)
    plt.plot(voltages[rpc][3:], good[rpc], label=labels[rpc]+" in-peak", marker = "o", color=colour[rpc])
    plt.plot(voltages[rpc][3:], bad[rpc], label=labels[rpc]+" off-peak", marker = "x", linestyle = "--", color=colour[rpc])

plt.title("Hit Counts vs. Voltage")
plt.ylabel("hit Count per Chunk")
plt.xlabel("Voltage / V")
plt.xlim(5400, 6000)
plt.ylim(0, 1200)
plt.legend()
plt.show()

labels = ["bottom triplet", "middle triplet", "singlet"]
plt.figure(figsize=(10, 6))  # Create a new figure
for rpc in range(3):
    print(rpc)
    plt.plot(voltages[rpc][3:], cluster[rpc], label=labels[rpc], marker = "o", color=colour[rpc])
    plt.errorbar(voltages[rpc][3:], cluster[rpc], yerr=error[rpc],  marker = "o", color=colour[rpc])
plt.title("Cluster Size vs. Voltage")
plt.ylabel("Average Cluster size")
plt.xlim(5400, 6000)
plt.xlabel("Voltage / V")
plt.legend()
plt.show()

labels = ["bottom triplet", "middle triplet", "singlet"]
plt.figure(figsize=(10, 6))  # Create a new figure
for rpc in range(3):
    print(rpc)
    plt.plot(voltages[rpc][3:], eff[rpc], label=labels[rpc], marker = "o", color=colour[rpc])
plt.title("Maximal Reconstruction Efficiency vs. Voltage")
plt.ylabel("Efficiency")
plt.xlim(5400, 6000)
plt.xlabel("Voltage / V")
plt.legend()
plt.show()


labels = ["bottom triplet", "middle triplet", "singlet"]
plt.figure(figsize=(10, 6))  # Create a new figure
for rpc in range(3):
    print(rpc)
    plt.plot(voltages[rpc], stds[rpc], label=labels[rpc], marker = "o", color=colour[rpc])
plt.title("Standard deviation of hit time vs. Voltage")
plt.ylabel("$\sigma$ / (25/32) ns")
plt.xlim(5000, 6000)
plt.xlabel("Voltage / V")
plt.legend()
plt.show()
"""
