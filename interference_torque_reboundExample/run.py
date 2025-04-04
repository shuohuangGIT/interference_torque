import subprocess
from multiprocessing import Pool
import itertools
import gc
import numpy as np

# C
C_EXECUTABLE = "./rebound"

def run_simulation(args):
    dynamic_tau_nAeB, dynamic_param, mass_ratio, phi_angle, tau_nA0 = args
    cmd = [C_EXECUTABLE, str(dynamic_tau_nAeB), str(dynamic_param), str(mass_ratio), str(phi_angle), str(tau_nA0)]
    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        print(f"complete: dynamic_tau_nAeB={dynamic_tau_nAeB}, dynamic_param={dynamic_param}, mass_ratio={mass_ratio}, phi_angle={phi_angle}, tau_nA0={tau_nA0}")
    except subprocess.CalledProcessError as e:
        print(C_EXECUTABLE, str(dynamic_tau_nAeB), str(dynamic_param), str(mass_ratio), str(phi_angle), str(tau_nA0))
        print(f"error: dynamic_tau_nAeB={dynamic_tau_nAeB}, dynamic_param={dynamic_param}, mass_ratio={mass_ratio}, phi_angle={phi_angle}, tau_nA0={tau_nA0} | error message: {e.stderr}")

if __name__ == "__main__":

    param_combinations = itertools.product(
        np.logspace(1, 3, 10),          # dynamic_tau_nAeB, default 0
        np.logspace(-2,0,10)+1,          # dynamic_param, default 0
        [1.],          # mass_ratio, default 1
        [-1],           # phi_angle, default 0
        [0.5]          # tau_nA/tao_o, default: np.linspace(1,5,5)
    )
    ncores = 200
    with Pool(processes=200) as pool:
        pool.map(run_simulation, param_combinations)

    pool.close()
    pool.join()
    gc.collect()

    print("all complete")
