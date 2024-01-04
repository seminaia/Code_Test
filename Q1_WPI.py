import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Load data
data = pd.read_excel('~/Downloads/data-1.xlsx') # Change to the correct directory

# Birch-Murnaghan equation
def BM(V, a, b, c, d):
    return a + b * V**(-2/3) + c * V**(-4/3) + d * V**(-6/3)

# Extract volume data
Vol = data['Volume(A^3)']


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

# Figure 1 
plt.figure(figsize=(8, 6))
plt.title("Birch-Murnaghan Fit to Energy Data")  

for composition, popt in zip(list(data)[1:], popt_list):
    plt.scatter(Vol, data[composition], label=composition)
    fit = BM(Vol, *popt)  # Use the same `volume` data for consistency
    plt.plot(Vol, fit, label=f'Birch Murnaghan Fit ({composition})')
    
plt.ylabel('Energy(eV/atom)')
plt.xlabel('Volume(A^3)')
plt.legend()

# Figure 2
plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.title("Minimum Energy(eV/atom) vs Volume (A^3)")  # Updated title

# Directly use the existing min_energies matrix
e0 = min_energies[:, 1]  # Extract minimum energies
V0= min_energies[:,0]
plt.scatter(V0, e0)

# Mixing Energy
# E_mix=(E_Hf+E_Zr-E_mixed)/2
E1=min_energies[0,1] # Equilibrium Energy of pure HfTiO4 Composition
E2=min_energies[4,1] # Equilibrium Energy of pure ZrTiO4 Composition 
E_mixed=min_energies[2,1] # Equilibrium Energy at 0.5 HfTiO4 and 0.5 ZrTiO4 Composition
E_mix= (E1+E2-E_mixed)/2
print(f"Mixing Energies'{E_mix}'")


plt.xlabel('Volume (A^3')
plt.ylabel('Minimum Energies (eV/atom)') #a grid for clarity
plt.show()
