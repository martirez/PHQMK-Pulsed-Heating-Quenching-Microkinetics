import cantera as ct
import numpy as np
import os
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from multiprocessing import Pool
import pandas as pd
import seaborn as sns

        
plt.rcParams['axes.linewidth'] = 3
plt.rc('xtick', labelsize=30)
plt.rc('ytick', labelsize=30)
plt.rc('axes', labelsize=30)
plt.rc('legend', fontsize=24)
plt.rcParams['lines.markersize'] = 5
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['xtick.major.size']=10
plt.rcParams['xtick.major.width']=2
plt.rcParams['xtick.minor.size']=10
plt.rcParams['xtick.minor.width']=10
plt.rcParams['xtick.minor.size']=10
plt.rcParams['xtick.minor.width']=2
plt.rcParams['ytick.major.size']=10
plt.rcParams['ytick.major.width']=2
plt.rcParams['legend.edgecolor']='k'
plt.rcParams['axes.unicode_minus']=False
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['xtick.major.pad'] = 8
plt.rcParams['ytick.major.pad'] = 8
plt.rcParams['legend.handletextpad']=0.4
plt.rcParams['legend.columnspacing']=0.5
plt.rcParams['legend.labelspacing']=0.3
plt.rcParams['legend.title_fontsize'] = 14
plt.rcParams["errorbar.capsize"] = 6

results = []

def run_PHQ(T_profile, output_filename, para):
    global results
    #######################################################################
    # Input Parameters
    #######################################################################
    volumetric_flow=1/1e6 #volumetric flow rate in m3/s 
    pressure=101325 #pressure Pa
    reactor_volume=1/1e6 #m^3
    specific_cat_area = 150000  # Catalyst specific surface area in m-1

    #import temperature profile
    time_arr = T_profile[:,0]
    temp_arr = T_profile[:,1]
    temp = temp_arr[0]
    temp_arr_in = np.array([temp]*(len(temp_arr)))

    #input file containing the surface reaction mechanism
    cti_file1 ='thermo1.yaml'
    cti_file2 = 'thermo2.yaml'
    #import the gas phase model
    gas1 = ct.Solution(cti_file1, 'gas')
    gas1.TPX = temp, pressure, 'H2:0.75, N2:0.25'
    gas2 = ct.Solution(cti_file2, 'gas')
    gas2.TPX = temp, pressure, 'H2:0.75, N2:0.25'
    bulk1 = ct.Solution(cti_file1, 'bulk')
    bulk1.TP = temp, pressure
    bulk2 = ct.Solution(cti_file2, 'bulk')
    bulk2.TP = temp, pressure
    # import the surface model
    surf1 = ct.Interface(cti_file1,'surface1', [gas1, bulk1])
    surf1.TP = temp, pressure
    surf1.coverages = {'FE(S3)': 0.667, 'N(S2)': 0.222, 'FE(S1)':0.111}
    surf2 = ct.Interface(cti_file2,'surface1', [gas2, bulk2])
    surf2.TP = temp, pressure
    surf2.coverages = {'FE(S3)': 0.889, 'N(S2)':0.111}
    # Loop through reactions and print their rate types
    for i, r in enumerate(surf1.reactions()):
        print(f"Reaction {i}: {r.equation}")
        print(f" → Rate type: {type(r.rate)}\n")
        
    ####################################################################
    rvol = reactor_volume

    ## catalyst area in one reactor
    cat_area = rvol*specific_cat_area
    mass_flow_rate =  volumetric_flow*gas1.mean_molecular_weight*pressure/ct.gas_constant/temp
    mdot_data =  (volumetric_flow * gas1.mean_molecular_weight * pressure) / (ct.gas_constant * temp_arr_in)
    mdot_func = interp1d(time_arr, mdot_data, kind='linear', fill_value='extrapolate')
    initial_molar_flux=mass_flow_rate/gas1.mean_molecular_weight*1000

    # Create well-stirred reactor and reservoirs (upstream and downstream from reactor)
    r1 = ct.IdealGasReactor(gas1, energy='off')
    r2 = ct.IdealGasReactor(gas2, energy='off')
    r1.volume = rvol
    r2.volume = rvol
    upstream1 = ct.Reservoir(gas1, name='upstream1')
    upstream2 = ct.Reservoir(gas2, name='upstream2')
    downstream1 = ct.Reservoir(gas1, name='downstream1')
    downstream2 = ct.Reservoir(gas2, name='downstream2')

    # Initialize catalyst surface in reactor
    rsurf1 = ct.ReactorSurface(surf1, r1, A=cat_area)
    rsurf2 = ct.ReactorSurface(surf2, r2, A=cat_area)
        
    # Add upstream (inlet) flow
    m1 = ct.MassFlowController(upstream1, r1, mdot=mdot_func)
    m2 = ct.MassFlowController(upstream2, r2, mdot=mdot_func)
    # Add downstream (outlet) flow
    v1 = ct.MassFlowController(r1, downstream1, mdot=mdot_func)
    v2 = ct.MassFlowController(r2, downstream2, mdot=mdot_func)
    
    # the simulation only contains one reactor    
    sim1 = ct.ReactorNet([r1])
    sim2 = ct.ReactorNet([r2])
    # The maximum number of error test failures permitted by the CVODES integrator in a single step. 
    sim1.max_err_test_fails = 100
    sim2.max_err_test_fails = 100

    # tolerances
    sim1.rtol = 1.0e-8
    sim1.atol = 1.0e-10
    sim2.rtol = 1.0e-8
    sim2.atol = 1.0e-10

    #set arrays
    mole_fracs1=np.zeros((len(time_arr),gas1.n_species))
    mole_fracs2=np.zeros((len(time_arr),gas2.n_species))
    coverages1=np.zeros((len(time_arr),surf1.n_species))
    coverages2=np.zeros((len(time_arr),surf2.n_species))
    weighted_results1 = np.zeros((len(time_arr), surf1.n_species))
    weighted_results2 = np.zeros((len(time_arr), surf2.n_species))
    molar_flux=np.zeros(len(time_arr))
    increments1=np.zeros(len(time_arr))
    increments2=np.zeros(len(time_arr))
    t_increments = np.zeros(len(time_arr))
    mass_flow = np.zeros(len(time_arr))
    mass_flow2 = np.zeros(len(time_arr))
    weights_array1 = np.zeros(len(time_arr))
    weights_array2 = np.zeros(len(time_arr))
    new_pressure = 101325
    
    #temperature pulsing
    for b in range(len(time_arr)): 
        if b == 0:
            gas1.TPX = temp, new_pressure, 'H2:0.75, N2:0.25'
            gas2.TPX = temp, new_pressure, 'H2:0.75, N2:0.25'
        gas1.TPX = temp, new_pressure, gas1.X
        gas2.TPX = temp, new_pressure, gas2.X
        surf1.TPX = temp, new_pressure, surf1.X
        surf2.TPX = temp,  new_pressure, surf2.X
        r1.syncState()
        r2.syncState()
        temp = temp_arr[b]
        gas1.TPX = temp, pressure, 'H2:0.75, N2:0.25'
        gas2.TPX = temp, pressure, 'H2:0.75, N2:0.25'
        upstream1.syncState()
        upstream2.syncState()
        sim1.advance(time_arr[b])
        sim2.advance(time_arr[b])
        increments1[b] = sim1.time   # time 
        increments2[b] = sim2.time   # time 
        coverages1[b,:]=surf1.X      # surface coverages
        coverages2[b,:]=surf2.X      # surface coverages
        mole_fracs1[b,:]=gas1.X      # gas mole fractions
        mole_fracs2[b,:]=gas2.X      # gas mole fractions
        t_increments[b]= temp      # temp at time step     
        if b == 0: 
           weights_array1[b] = 0.864
           weights_array2[b] = (1-weights_array1[b])
           weighted_results1[b, :] = coverages1[b, :] * weights_array1[b]
           weighted_results2[b, :] = coverages2[b, :] * weights_array2[b]
           surf1.X = weighted_results1[b, :] + weighted_results2[b, :]
           surf2.X = weighted_results1[b, :] + weighted_results2[b, :]
           gas1.X = mole_fracs1[b,:]* weights_array1[b] + mole_fracs2[b,:]* weights_array2[b]
           gas2.X = mole_fracs1[b,:]* weights_array1[b] + mole_fracs2[b,:]* weights_array2[b]
        if b > 0: 
           # Species to exclude
           species_duplicate = ["N2(S2)", "H2(S1)", "NH3(S1)"]  
           # Get species names
           species_names1 = surf1.species_names
           species_names2 = surf1.species_names
           # Find indices of species to keep
           keep_indices1 = [i for i, species in enumerate(species_names1) if species not in species_duplicate]
           keep_indices2 = [i for i, species in enumerate(species_names1) if species not in species_duplicate]
           # Filter species names
           species_names1_filtered = [species_names1[i] for i in keep_indices1]
           species_names2_filtered = [species_names2[i] for i in keep_indices2]
           coverages1_new = weighted_results1[b-1, keep_indices1] + weighted_results2[b-1, keep_indices2]
           coverages1_r = coverages1_new / np.sum(coverages1_new)
           renormalized_values1 = dict(zip(species_names1_filtered, coverages1_r))
           desired_species = ["FE(S3)", "FE(S1)"] 
           surf_tot = 1-(renormalized_values1.get("FE(S3)", 0) + renormalized_values1.get("FE(S1)", 0))
           weights_array1[b] = 1/(1+np.exp(-(surf_tot-(3/18))*100/3))
           weights_array2[b] = (1-weights_array1[b])
           weighted_results1[b, :] = coverages1[b, :] * weights_array1[b]
           weighted_results2[b, :] = coverages2[b, :] * weights_array2[b]
           surf1.X = weighted_results1[b, :] + weighted_results2[b, :]
           surf2.X = weighted_results1[b, :] + weighted_results2[b, :]
           gas1.X = mole_fracs1[b,:]* (weights_array1[b]) + mole_fracs2[b,:]* (weights_array2[b])
           gas2.X = mole_fracs1[b,:]* (weights_array1[b]) + mole_fracs2[b,:]* (weights_array2[b])
        if b == 0: 
            mass1 = r1.mass
            mass2 = r2.mass 
        new_pressure_1 = mass1*ct.gas_constant*temp/gas1.mean_molecular_weight/reactor_volume  
        new_pressure_2 = mass2*ct.gas_constant*temp/gas2.mean_molecular_weight/reactor_volume  
        new_pressure = new_pressure_1*weights_array1[b] + new_pressure_2*weights_array2[b]

    species_duplicate = ["N2(S2)", "H2(S1)", "NH3(S1)"]  # Modify based on what you want to exclude
    species_names1 = surf1.species_names
    species_names2 = surf2.species_names
    keep_indices1 = [i for i, species in enumerate(species_names1) if species not in species_duplicate]
    keep_indices2 = [i for i, species in enumerate(species_names2) if species not in species_duplicate]
    species_names1_filtered = [species_names1[i] for i in keep_indices1]
    species_names2_filtered = [species_names2[i] for i in keep_indices2]
    coverages1_new = coverages1[:, keep_indices1]
    coverages2_new = coverages2[:, keep_indices2]
    weighted_results1_new = weighted_results1[:, keep_indices1]
    weighted_results2_new = weighted_results2[:, keep_indices2]
    coverages1_r = coverages1_new / coverages1_new.sum(axis=1, keepdims=True)
    coverages2_r = coverages2_new / coverages2_new.sum(axis=1, keepdims=True)
    weights_array1_1 = np.tile(weights_array1[:, np.newaxis], (1, 11))
    weights_array2_1 = np.tile(weights_array2[:, np.newaxis], (1, 11))
    weights_array1_2 = np.tile(weights_array1[:, np.newaxis], (1, 3))
    weights_array2_2 = np.tile(weights_array2[:, np.newaxis], (1, 3))
    weighted_results1_r = coverages1_new / coverages1_new.sum(axis=1, keepdims=True)*weights_array1_1
    weighted_results2_r = coverages2_new / coverages2_new.sum(axis=1, keepdims=True)*weights_array2_1
    mole_fracs1_weighted=mole_fracs1*weights_array1_2
    mole_fracs2_weighted=mole_fracs2*weights_array2_2
    
    
    desired_species = ["N(S2)", "N(S1)", "NH(S1)", "FE(S3)", "FE(S1)"]  # Modify based on the species of interest

    # Get indices for desired species
    desired_indices1 = [species_names1_filtered.index(species) for species in desired_species if species in species_names1_filtered]
    desired_indices2 = [species_names2_filtered.index(species) for species in desired_species if species in species_names2_filtered]
    
    # Initialize filtered storage arrays
    filtered_coverages1 = np.zeros((len(time_arr), len(desired_species)))
    filtered_coverages2 = np.zeros((len(time_arr), len(desired_species)))
    filtered_wcoverages1 = np.zeros((len(time_arr), len(desired_species)))
    filtered_wcoverages2 = np.zeros((len(time_arr), len(desired_species)))

    # Loop through time steps and extract coverages
    for b in range(len(time_arr)):
        filtered_coverages1[b, :] = coverages1_r[b, desired_indices1]
        filtered_coverages2[b, :] = coverages2_r[b, desired_indices2]
        filtered_wcoverages1[b, :] = weighted_results1_r[b, desired_indices1]
        filtered_wcoverages2[b, :] = weighted_results2_r[b, desired_indices2]
       
    T = t_increments
    sorted_coverages = filtered_wcoverages1[(T > 200.15) & (T < 1473.15)] + filtered_wcoverages2[(T > 200.15)& (T < 1473.15)] 
    sorted_moles = mole_fracs1_weighted[(T > 200.15) & (T < 1473.15)] + mole_fracs2_weighted[(T > 200.15)& (T < 1473.15)] 
    sorted_time  = increments1[(T > 200.15) & (T < 1473.15)]     
    
    data=np.c_[sorted_time, t_increments, sorted_moles, sorted_coverages, weights_array1, weights_array2]
    print('writing results into {}'.format(output_filename))
    np.savetxt(output_filename, data, header = 'x_s, T, H2, N2, NH3, N(S2), N(S1), NH(S1), FE(S3), FE(S1), w1, w2', delimiter=",")
    x_s = data[:, 0]  # time (x_s)
    N_S2 = data[:, 5]  # 'N(S2)' column
    N_S1 = data[:, 6]
    NH_S1 = data[:, 7]  # 'NH(S1)' column
    FE_S3 = data[:, 8]  # 'FE(S3)' column
    N_tot = N_S2 + N_S1
    FE_S1 = data[:, 9]  # 'FE(S1)' column
    FE_tot = FE_S3 + FE_S1
    H2 = data[:, 2]
    N2 = data[:, 3]
    NH3 = data[:, 4]
    T = data[:, 1] 
    
    sort = np.argsort(x_s) #sorted indices

    # coverages
    plt.figure(figsize=(10, 10))
    plt.plot(x_s[sort], N_tot[sort], linewidth=5.0,  color = 'rebeccapurple', label='$\mathrm{*N}$')
    plt.plot(x_s[sort], NH_S1[sort], linewidth=5.0, color = 'lightseagreen', label='$\mathrm{*NH}$')
    plt.plot(x_s[sort], FE_tot[sort], linewidth=5.0, color = 'mediumvioletred', label='$\mathrm{*}$')
    plt.xlabel("time [s]")
    plt.ylabel("coverage [ML]")
    plt.legend()
    plt.xlim([10,15])
    plt.ylim([0,1])
    plt.grid(False)

    plt.savefig(f"{output_filename}_cov.png", dpi=300, bbox_inches="tight")

    # gas species 
    plt.figure(figsize=(10, 10))

    #plt.plot(x_s[sort], H2[sort], linewidth=5.0,  color = 'red', label='$\mathrm{H_2}$')
    #plt.plot(x_s[sort], N2[sort], linewidth=5.0, color = 'purple', label='$\mathrm{N_2}$')
    plt.plot(x_s[sort], NH3[sort]*1e6, linewidth=5.0, color = 'seagreen', label='$\mathrm{NH_3}$')
    NH3_sort = NH3[sort]
    x_s = x_s[sort]

    plt.xlabel("time [s]")
    plt.ylabel(r'mole fraction NH$_3$ [ppm]')
    plt.legend()
    plt.xlim([10,15])
    plt.ylim([0,1000])
    plt.grid(False)
    plt.savefig(f"{output_filename}_gas.png", dpi=300, bbox_inches="tight")
    
    #temperature profile
    plt.figure(figsize=(10, 10))
    plt.plot(x_s, T, linewidth=5.0, color = 'mediumvioletred')
    plt.xlabel('$\mathrm{time\ [s]}$')
    plt.ylabel('$\mathrm{temperature\ [K]}$')
    plt.xlim([10,15])
    plt.ylim([350,1450])
    plt.title('$\mathbf{Temperature\ Pulsing}$',fontsize=30)
    plt.savefig(f"{output_filename}_temp.png", transparent=False, bbox_inches='tight')


    #average_ammonia
    window_1 = (x_s >= 14.5) & (x_s <= 15.5)
    window_2 = (x_s >= 22.5) & (x_s <= 23.5)

    # Get time and temperature within those windows
    x_1 = x_s[window_1]
    x_2 = x_s[window_2]
    T_1 = T[window_1]
    T_2 = T[window_2]
    NH3_1 = NH3_sort[window_1]
    NH3_2 = NH3_sort[window_2]
    idx_min_1 = np.argmin(T_1)
    idx_min_2 = np.argmin(T_2)
    t_min_1 = x_1[idx_min_1]
    t_min_2 = x_2[idx_min_2]
    mask_between_mins = (x_s >= min(t_min_1, t_min_2)) & (x_s <= max(t_min_1, t_min_2))
    NH3_between = NH3_sort[mask_between_mins]
    T_between = T[mask_between_mins]

    # Compute averages
    avg_NH3_between = np.mean(NH3_between)
    avg_T_between = np.mean(T_between)

    # Store results
    results.append({
        'Th': para[0],
        'Tc': para[1],
        'avg_NH3_between_mins': avg_NH3_between,
        'avg_T_between_mins': avg_T_between
    })

    print("Just added to results:", results[-1])

    return {
        'Th': para[0],
        'Tc': para[1],
        'avg_NH3_between_mins': avg_NH3_between,
        'avg_T_between_mins': avg_T_between
    }
    
def cal_T_profile(T_h, T_c, para_h, para_rc, para_cc, dt=0.0001, t_total=25, T_th=1100):
    #generate pulse heating temperature profile T_profile=[t, T]
    #T_h: heating(maximum) temperature (K)
    #T_c: cooling(original) temperature (K)
    #para_h: parameter for heating zone
    #para_rc: parameters for radiation cooling zone
    #para_cc: parameters for convection cooling zone
    #T_th: temperature thresold, where radiation cooling transform into convection cooling
    #dt: time step (s)
    #t_total: total duration

    t = np.arange(0, t_total, dt)
    T = np.zeros(len(t))
    cyc = 1  # duration of one pulse in seconds
    pulse_len = int(cyc / dt)
    num_cycles = len(t) // pulse_len

    current_T = T_c  # initial temperature

    full_profile = []

    for cycle in range(num_cycles):
        pulse_T = np.zeros(pulse_len)
        pulse_T[0] = current_T
        i = 0

        # Heating
        while i + 1 < pulse_len and pulse_T[i] < T_h:
            pulse_T[i + 1] = pulse_T[i] + dt * para_h
            i += 1

        # Radiative cooling
        while i + 1 < pulse_len and pulse_T[i] > T_th:
            pulse_T[i + 1] = pulse_T[i] + dt * (-para_rc[0] * pulse_T[i]**4 + para_rc[1])
            i += 1

        # Convective cooling
        while i + 1 < pulse_len:
            pulse_T[i + 1] = pulse_T[i] + dt * (para_cc[0] * pulse_T[i] + para_cc[1])
            i += 1

        # Extend convective cooling if needed
        while len(pulse_T) < pulse_len:
            next_T = pulse_T[-1] + dt * (para_cc[0] * pulse_T[-1] + para_cc[1])
            pulse_T = np.append(pulse_T, next_T)

        full_profile.extend(pulse_T)
        current_T = pulse_T[-1] # update starting temp for next cycle

    # Trim to match array length
    t = np.array(t[:len(full_profile)])
    T_profile = np.vstack((t, full_profile)).T

    return T_profile

from scipy.optimize import minimize_scalar

def cal_pulse_error(fre, para, para_h_ori, para_rc_ori, para_cc_ori, target_cycle_duration, T_th=1100, dt=0.0001):
    Tc = para[1]
    
    # Generate temperature profile with current frequency factor
    T_profile_trial = cal_T_profile(
        para[0], Tc,
        para_h_ori * fre,
        para_rc_ori(Tc) * fre,
        para_cc_ori(Tc) * fre,
        dt=dt, t_total=25, T_th=T_th
    )
    
    # Find pulse cycle length: when it returns near Tc
    tol = 0.05
    matches = np.where(np.abs(T_profile_trial[:, 1] - Tc) < tol)[0]
    if len(matches) >= 2:
        pulse_time = T_profile_trial[matches[1], 0]
        return abs(pulse_time - target_cycle_duration)  # Goal is to minimize difference from 1 second
    else:
        return np.inf  # invalid cycle

def single_task(para):
    durations = [1.0]
    Tc = para[1]
    para_h_ori = 3105
    para_rc_ori = lambda Tc: np.array([5.039e-9, 5274.4+(1.912*Tc/1.0657)])
    para_cc_ori = lambda Tc: np.array([-1.912, 1.912*Tc/1.0657])
    T_profile_final = cal_T_profile(para[0],para[1],para_h_ori*para[2],para_rc_ori(para[1])*para[3],para_cc_ori(para[1])*para[3], dt=0.0001, t_total=25)
    current_dir = os.getcwd()
    name = "Th={}_Tc={}_hr={}_cr={}.csv".format(para[0],para[1],para[2],para[3])
    filename = os.path.join(current_dir,name)
    metrics = run_PHQ(T_profile_final,filename, para)
    return metrics

if __name__ == '__main__':
    para_h_ori = 3105
    para_rc_ori = lambda Tc: np.array([5.039e-9, 5274.4+(1.912*Tc/1.0657)])
    para_cc_ori = lambda Tc: np.array([-1.912, 1.912*Tc/1.0657])

    base_values = [1332]
    T_h_arr = np.array(base_values)
    base_values = [925]
    T_c_arr = np.array(base_values)
    fre_factor_arr = (np.linspace(1, 2, 1))
    combinations = [(Th, Tc, fre, fre) for Th,Tc in zip(T_h_arr,T_c_arr) for fre in fre_factor_arr]
    print (combinations)
    with Pool(processes=8) as pool:
        results = pool.map(single_task, combinations)
        
print("Length of results:", len(results))
for i, entry in enumerate(results[:5]):
    print(f"Entry {i}: {entry}")
df = pd.DataFrame(results)
print("Columns in df:", df.columns)
print("Sample rows:\n", df.head())
df.to_csv("pulse_metrics.csv", index=False)
pivot = df.pivot_table(index='Th', columns='Tc', values='avg_NH3_between_mins')
plt.figure(figsize=(10, 8))
ax = sns.heatmap(pivot, cmap="YlGnBu", fmt=".4f", vmin = 0, vmax = 0.001, cbar_kws={'label': 'mole fraction NH$_3$'})
ax.tick_params(axis='both', direction='out')
plt.gca().invert_yaxis()
plt.xlabel("Min Temperature [K]")
plt.ylabel("Max Temperature [K]")
plt.tight_layout()
plt.savefig("heatmap_NH3_vs_Textremes.png", dpi=300)

