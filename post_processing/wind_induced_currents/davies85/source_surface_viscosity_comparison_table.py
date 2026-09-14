import numpy as np

S = 2.
G = 0.3
Cd = 2.5*10**-3

k0 = 0.4
f = 1.2 * 10**-4

g = 9.81

rhoA = 1.225
rho = 1025

w = np.array([5., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30.]) 

tau = Cd * rhoA * w**2
u = np.sqrt(tau/rho)

h = k0*u/f

nuA = k0*S*u*G*w**2/g
nuB = (0.1825 * 10**-3) * w**2.5

wp = w * 1.95**0.17

H = 0.0213 * wp**2
T = 0.5635 * wp

nuC = (0.028 * H**2) / T

output.RowData.append(w, f'Wind Speed, W')
output.RowData.append(tau, f'Wind Stress')
output.RowData.append(u, f'Frictional velocity')
output.RowData.append(h, f'Depth of Frictional Influence')
output.RowData.append(nuA, f'nuA')
output.RowData.append(nuB, f'nuB')
output.RowData.append(nuC, f'nuC')
