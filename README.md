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
replace Tmin (T_c_arr) and Tmax (T_h_arr) 
  
Execute:  
conda activate ct-env  
#for isothermal stead state  
python Fe-NH3_iso.py  
or  
#for pulsed heating and quenching mode  
python Fe-NH3_phq.py  

