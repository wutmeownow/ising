import random
import math
import time
import os

# Parameters
NX = 64
NY = 64
ntherm = 0
VisualDisplay = 1
SleepTime = 100000  # in microseconds

# Initialize the lattice with boundary spins
spin = [[0 for _ in range(NY + 2)] for _ in range(NX + 2)]

def update_spin(s, env):
    """Do a metropolis update on a spin s whose environment is env"""
    spin = s[0]
    newspin = 1 if random.random() < 0.5 else -1
    DeltaBetaE = -(newspin - spin) * env
    if DeltaBetaE <= 0 or random.random() < math.exp(-DeltaBetaE):
        s[0] = newspin

def sweep(beta, h):
    """Sweep through all lattice sites"""
    for nx in range(1, NX + 1):
        for ny in range(1, NY + 1):
            environment = (beta * (spin[nx][ny-1] + spin[nx][ny+1] + 
                                 spin[nx-1][ny] + spin[nx+1][ny]) + h)
            # Using list to simulate pointer
            spin_ref = [spin[nx][ny]]
            update_spin(spin_ref, environment)
            spin[nx][ny] = spin_ref[0]

def initialize_hot():
    """Initialize lattice with random spins"""
    print("Initializing system")
    for nx in range(NX + 2):
        for ny in range(NY + 2):
            if nx == 0 or nx == NX + 1 or ny == 0 or ny == NY + 1:
                spin[nx][ny] = 0
            else:
                spin[nx][ny] = 1 if random.random() < 0.5 else -1

def magnetization():
    """Calculate average magnetization"""
    nmag = sum(spin[nx][ny] for nx in range(1, NX + 1) 
              for ny in range(1, NY + 1))
    return nmag / (NX * NY)

def display_lattice(sweep_number):
    """Display the lattice configuration"""
    if SleepTime > 0:
        os.system('cls' if os.name == 'nt' else 'clear')
    
    for nx in range(1, NX + 1):
        for ny in range(1, NY + 1):
            print('X' if spin[nx][ny] == 1 else '-', end='')
        print()
    
    if sweep_number >= 0:
        print(f"sweep {sweep_number}:   magnetization <sigma> = {magnetization():.6f}")
    
    if SleepTime > 0:
        time.sleep(SleepTime / 1_000_000)  # Convert microseconds to seconds
    else:
        print()

def main():
    print(f"Program generates a thermal ensemble a 2D Ising model of "
          f"{NX}x{NY} spins with free boundary conditions.\n")
    
    random.seed(int(time.time()))
    
    nsweep = int(input("Enter total number of configurations generated:\n"))
    h = float(input("Enter value of magnetic field parameter h:\n"))
    beta = float(input("Enter temperature parameter beta (= 1/kT):\n"))
    
    initialize_hot()
    
    # Thermalization sweeps
    print(f"Thermalizing system, {ntherm} sweeps")
    for n in range(ntherm):
        sweep(beta, h)
        if VisualDisplay and n % (ntherm//10 if ntherm > 0 else 1) == 0:
            display_lattice(-1)
            print(f"Thermalization sweep {n}")
            time.sleep(1)
    
    # Main sweeps
    nmag = ntotal = 0
    for n in range(nsweep):
        if VisualDisplay:
            display_lattice(n)
        sweep(beta, h)
        for nx in range(1, NX + 1):
            for ny in range(1, NY + 1):
                nmag += spin[nx][ny]
                ntotal += 1
    
    if VisualDisplay:
        display_lattice(nsweep)
    print(f"Average Magnetization: <s> = {nmag/ntotal:.6f}")

if __name__ == "__main__":
    main()
    
