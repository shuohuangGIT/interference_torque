# interference_torque
The implement of interference torque (Yang&amp;Li 2024) based on rebound and reboundx. If you use the code, please cite xx.  

Put the whole file in the 'reboundx/examples/' directory. The way to run the default simulations:  
make clean  
make # read and compile problem.c where interference torque is implemented  
python3 run.py # run the simulations  

run.py is the python interface, where you can change the input parameters and also the number of CPU processors (ncores) that will be used. 
