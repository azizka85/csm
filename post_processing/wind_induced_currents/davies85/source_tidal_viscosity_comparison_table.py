import numpy as np

k = 2. * 10**-5
sigma = 1.2 * 10**-4

ut = np.array([0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0])

mum = k * ut**2 / sigma

output.RowData.append(ut, f'Tidal Velocity, UT')
output.RowData.append(mum, f'Tidal Viscosity, MuT')
