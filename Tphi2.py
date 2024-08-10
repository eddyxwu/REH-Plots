import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from scipy import special
from scipy.optimize import curve_fit

temp = np.array([13.79, 15.2, 36.45, 47.94, 50.77, 57.21, 69.97, 73.12, 88.12, 97.82, 114.85, 119.45, 138.25, 140.28, 156.03, 168.69, 184.98]) * 1e-3
temp_mod = np.array([15.2, 36.45, 50.77, 69.97, 73.12, 88.12, 97.82, 114.85, 119.45, 138.25, 140.28, 156.03, 168.69, 184.98]) * 1e-3

Q2_Tphi = np.array([99.97, 89.29, 90.09, 88.91, 70.83, 73.52, 59.72, 49.41, 36.18, 24.31, 15.18, 14.85, 9.34, 10.92, 7.96, 6.63, 10.96]) * 1e-6
Q3_Tphi = np.array([81.70, 83.06, 62.02, 46.46, 47.99, 30.14, 25.10, 16.47, 16.61, 8.92, 11.27, 8.44, 6.84, 7.21]) * 1e-6
Q4_Tphi = np.array([75.54, 67.51, 60.49, 47.69, 39.99, 28.76, 25.27, 19.22, 15.88, 11.31, 11.51, 8.60, 8.25, 5.91]) * 1e-6

plt.plot(temp, Q2_Tphi, marker = 'o', color = 'b', label='Q2')
plt.plot(temp_mod, Q3_Tphi, marker = 'o', color='red', label='Q3')
plt.plot(temp_mod, Q4_Tphi, marker = 'o', color='orange', label='Q4')

plt.title('T_phi as temperature is raised')
plt.xlabel('Temperature (K)')
plt.ylabel("Time (sec)")

plt.legend()
plt.show()
