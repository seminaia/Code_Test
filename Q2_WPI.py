import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Load data
data = pd.read_csv('~/Downloads/eos.csv')

# Birch-Murnaghan equation
def BM(V, a, b, c, d):
    return a + b * V**(-2/3) + c * V**(-4/3) + d * V**(-6/3)

# Extract volume data
Vol = data['Volume (A^3/atom)']

# Create a matrix to store minimum energies (outside of the loop)
min_energies = np.zeros((len(list(data)[1:]), 2))  # 5 rows for 5 compositions

# Perform fitting for each composition
popt_list = []

for i, composition in enumerate(list(data)[1:]):
    energy = data[composition]

    # Fit the Birch-Murnaghan equation
    popt, pcov = curve_fit(BM, Vol, energy)
    popt_list.append(popt)

    # Print fitting parameters
    print(f"Composition '{composition}' parameters:")
    print(f"\t a = {popt[0]:.4f}, b = {popt[1]:.4f}, c = {popt[2]:.4f}, d = {popt[3]:.4f}")

    # Generate the fitted curve
    fit = BM(Vol, *popt)

    # Find the index of the minimum energy point
    min_energy_idx = np.argmin(fit)

    # Extract the corresponding volume and energy
    min_volume = Vol[min_energy_idx]
    min_energy = fit[min_energy_idx]

    # Store volume and energy in the matrix
    min_energies[i, 0] = min_volume
    min_energies[i, 1] = min_energy

    # Errors
    residuals = energy - fit
    rmse = np.sqrt(np.mean(residuals**2))
    mae = np.mean(np.abs(residuals))

    # Output fitting errors
    print(f"Fitting errors for composition '{composition}':")
    print(f"\tRMSE: {rmse:.4f} eV/atom")
    print(f"\tMAE: {mae:.4f} eV/atom")
    print(f"\tE0: {min_energies[i, 1]:.4f} eV/atom at volume {min_energies[i, 0]:.4f} A^3")

 
plt.figure(figsize=(8, 6))
plt.title("Birch-Murnaghan Fit to Energy Data")

for composition, popt in zip(list(data)[1:], popt_list):
    plt.scatter(Vol, data[composition], label=composition)
    fit = BM(Vol, *popt)
    plt.plot(Vol, fit, label=f'Birch Murnaghan Fit ({composition})')
    
plt.ylabel('Energy(eV/atom)')
plt.xlabel('Volume(A^3/atom)')

# Directly use the existing min_energies matrix
e0 = min_energies[:, 1]  # Extract minimum energies
V0= min_energies[:,0]
plt.scatter(V0, e0)
plt.legend()
plt.show()
