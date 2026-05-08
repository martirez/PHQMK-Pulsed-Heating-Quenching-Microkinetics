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
#for pulsed heating and quenching mode  
python Fe-NH3_phq.py  
#for pulsed heating and quenching mode scaled rates  
python Fe-NH3_phq_scaled.py   

**Edit specifications in Fe-NH3_\<type\>.py:**  
replace Tmin (base_values, T_c_arr) and Tmax (base_values, T_h_arr) in K for PHQ    
replace temperature(s) (T_list) in K for isothermal      
replace pulse cycle length (cyc) in s for PHQ   
replace averaging window: lower bound (window_1) and upper bound (window_2) for all types     
replace simulation length (t_total) in s for all types 
replace cycle length to constrain Tmin (durations) in s for PHQ scaled rates
  
**Please cite this repository and the following publications:**  
**(1)** S. Kurdziel and J. M. P. Martirez, "PHQMK-Pulsed-Heating-Quenching-Microkinetics", (2026) https://github.com/martirez/PHQMK-Pulsed-Heating-Quenching-Microkinetics  
**(2)** S. Kurdziel, Z. Sun, Z. Chang, F. Liu, L. Hu, Y. Ju, E. A. Carter, and J. M. P. Martirez, "First-Principles Microkinetic Model for Ammonia Synthesis on Fe(110) under Non-Equilibrium Pulsed Heating and Quenching Operation", XXX, (2026)  
**(3)** J. M. P. Martirez, S. Kurdziel, and E. A. Carter, "First-Principles Insights into the Thermocatalytic Cracking of Ammonia-Hydrogen Blends on Fe(110). 2. Kinetics," Journal of Physical Chemistry C, 129, 6697-6718 (2025)  
**(4)** J. M. P. Martirez and E. A. Carter, First-Principles Insights into the Thermocatalytic Cracking of Ammonia-Hydrogen Blends on Fe(110). 1. Thermodynamics, Journal of Physical Chemistry C, 126, 19733 (2022)    
 
