import numpy as np

C1 = inputs[0].PointData['C']
C2 = inputs[1].PointData['C']

diff = np.abs(C1 - C2)

output.PointData.append(diff, f'diff')
