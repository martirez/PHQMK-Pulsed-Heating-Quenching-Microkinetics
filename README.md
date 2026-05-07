# PHQMK-Pulsed-Heating-Quenching-Microkinetics

**Requirements:**  
Cantera 3.1.0 (released December 2024)  
compiled using conda create --name ct-env --channel conda-forge cantera ipython matplotlib jupyter  
Tested only for 3.1.0, currently incompatible with 3.2.0 (released November 2025)
  
**Dependencies:**  
Anaconda 2023.09  
NumPy  
SciPy  
Matplotlib  
Pandas  
Seaborn  
  
**Input files:**   
thermo1.yaml (catalytic cycle A energetics)  
thermo2.yaml (catalytic cycle B energetics)  
gibbs_rate.py (external reaction type subroutine: surface-cubic-Ea)      
  
**Execute:**  
conda activate ct-env  
#for isothermal stead state  
python Fe-NH3_iso.py  
or  
#for pulsed heating and quenching mode  
python Fe-NH3_phq.py  
#for pulsed heating and quenching mode scaled-rates  
python Fe-NH3_phq_scaled.py   

**Edit specifications:**  
replace Tmin base_values (T_c_arr) and Tmax base_values (T_h_arr) in K for PHQ  
replace Ts (T_list) in K for isothermal    
replace pulse cycle length (cyc) in s  
replace averaging window: lower bound (window_1) and upper bound (window_2)   
replace simulation length (t_total) in s  
  
**Please cite this repository and the following publications:**  
(1) Sophia Kurdziel, Zijian Sun, Ziqiao Chang, Fangyuan Liu, Liangbing Hu, Yiguang Ju, Emily A. Carter,and John Mark P. Martirez, "First-Principles Microkinetic Model for Ammonia Synthesis on Fe(110) under Non-Equilibrium Pulsed Heating and Quenching Operation", XXX, (2026)  
(2) J. M. P. Martirez, S. Kurdziel, and E. A. Carter, "First-Principles Insights into the Thermocatalytic Cracking of Ammonia-Hydrogen Blends on Fe(110). 2. Kinetics," Journal of Physical Chemistry C, 129, 6697-6718 (2025)  
(3) J. M. P. Martirez and E. A. Carter, First-Principles Insights into the Thermocatalytic Cracking of Ammonia-Hydrogen Blends on Fe(110). 1. Thermodynamics, Journal of Physical Chemistry C, 126, 19733 (2022)    
 
