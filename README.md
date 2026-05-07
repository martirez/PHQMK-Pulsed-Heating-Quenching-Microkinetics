# PHQ-MK-Pulsed-Heating-Quenching-Microkinetics

Requirement:  
Cantera 3.1.0  
compiled using conda create --name ct-env --channel conda-forge cantera ipython matplotlib jupyter  
Tested only for 3.1.0, currently incompatible with 3.2.0 (released November 2025)
  
Dependencies:   
Anaconda 2023.09  
NumPy  
SciPy  
Matplotlib  
Pandas  
Seaborn  
  
Input files:   
thermo1.yaml (catalytic cycle A energetics)  
thermo2.yaml (catalytic cycle B energetics)  
gibbs_rate.py (external reaction type subroutine: surface-cubic-Ea)      
  
Edit specifications:  
replace Tmin base_values (T_c_arr) and Tmax base_values (T_h_arr) in K for PHQ  
replace T (T_list) in K for isothermal    
replace pulse cycle length (cyc) in s  
replace averaging window: lower bound (window_1) and upper bound (window_2)  
  
Execute:  
conda activate ct-env  
#for isothermal stead state  
python Fe-NH3_iso.py  
or  
#for pulsed heating and quenching mode  
python Fe-NH3_phq.py  
#for pulsed heating and quenching mode scaled-rates  
python Fe-NH3_phq_scaled.py   
  
Reference (please cite):  
Sophia Kurdziel, Zijian Sun, Ziqiao Chang, Fangyuan Liu, Liangbing Hu, Yiguang Ju, Emily A. Carter,and John Mark P. Martirez, "First-Principles Microkinetic Model for Ammonia Synthesis on Fe(110) under Non-Equilibrium Pulsed Heating and Quenching Operation", XXX, 2026  
 
