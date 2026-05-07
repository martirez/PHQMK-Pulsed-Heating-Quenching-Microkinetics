# PHQ-MK-Pulsed-Heating-Quenching-Microkinetics

Dependencies:   
Anaconda 2023.09  
Cantera 3.1.0  
NumPy  
SciPy  
Matplotlib  
Pandas  
Seaborn  
  
files: thermo1.yaml, thermo2.yaml, and gibbs_rate.py  
  
To use:  
replace Tmin base_values (T_c_arr) and Tmax base_values (T_h_arr) in K   
replace pulse cycle length (s): cyc  
  
Execute:  
conda activate ct-env  
#for isothermal stead state  
python Fe-NH3_iso.py  
or  
#for pulsed heating and quenching mode  
python Fe-NH3_phq.py  
#for pulsed heating and quenching mode scaled-rates  
python Fe-NH3_phq_scaled.py   
