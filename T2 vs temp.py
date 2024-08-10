import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# raw data
Q2_T2 = np.array([70.47, 63.14, 58.04, 60.73, 53.45, 43.64, 45.32, 35.77, 29.38, 20.93, 13.86, 13.61, 8.574, 10.12, 6.946, 5.864, 7.16]) * 1e-6
Q3_T2 = np.array([56.43, 50.33, 42.23, 36.91, 41.01, 24.03, 22.13, 14.93, 14.62, 8.404, 10.53, 7.709, 6.003, 4.743]) * 1e-6
Q4_T2 = np.array([48.96, 51.6, 46.14, 37.93, 32.72, 24.38, 20.98, 16.7, 14.79, 10.26, 10.57, 7.702, 7.041, 3.698]) * 1e-6

temp = np.array([13.79, 15.2, 36.45, 47.94, 50.77, 57.21, 69.97, 73.12, 88.12, 97.82, 114.85, 119.45, 138.25, 140.28, 156.03, 168.69, 184.98]) * 1e-3
temp_mod = np.array([15.2, 36.45, 50.77, 69.97, 73.12, 88.12, 97.82, 114.85, 119.45, 138.25, 140.28, 156.03, 168.69, 184.98]) * 1e-3

T2_avg = np.array([70.47, 56.17666667, 53.32333333, 60.73, 47.27333333, 43.64, 40.05333333, 36.5, 25.93, 21.34666667, 15.16333333, 14.34, 9.079333333, 10.40666667, 7.452333333, 6.302666667, 5.200333333])

# regular plot
plt.subplots()
plt.plot(temp, Q2_T2, marker = 'o', color = 'blue', label='Q2')
plt.plot(temp_mod, Q3_T2, marker = 'o', color='red', label='Q3')
plt.plot(temp_mod, Q4_T2, marker = 'o', color='orange', label='Q4')

plt.title('T2 as temperature is raised')
plt.xlabel('Temperature (K)')
plt.ylabel("T2 (sec)")

plt.legend()
plt.show()


# loglog plot
plt.subplots()
plt.loglog(temp, Q2_T2, marker = 'o', label='Q2')
plt.loglog(temp_mod, Q3_T2, marker = 'o', color='red', label='Q3')
plt.loglog(temp_mod, Q4_T2, marker = 'o', color='orange', label='Q4')

plt.title('T2 as temperature is raised')
plt.xlabel('Temperature (mK)')
plt.ylabel("T2 (sec)")

plt.legend()
plt.show()

# model fits
temp_fit = np.array([88.12, 97.82, 114.85, 119.45, 138.25, 140.28, 156.03, 168.69, 184.98])
T2_avg_fit = np.array([25.93, 21.34666667, 15.16333333, 14.34, 9.079333333, 10.40666667, 7.452333333, 6.302666667, 5.200333333])

# def func(temp_fit, a, b, c):
#     return a * np.exp(temp_fit / b) + c

# def func(temp_fit, a, b, c):
#     return a * temp_fit ** b + c
#
# popt, pcov = curve_fit(func, temp_fit, T2_avg_fit, p0 = [260, -0.16, -106], maxfev = 10000) # bounds=(0,10))
# popt
# a_opt, b_opt, c_opt = popt
# x_model = np.linspace(min(temp_fit), max(temp_fit)*5, 100)
# y_model = func(x_model, a_opt, b_opt, c_opt)
# plt.plot(x_model, y_model)
#
#
# residuals = T2_avg_fit - func(temp_fit, *popt)
# ss_res = np.sum(residuals**2)
# print('Sum of squares of the residuals: ', ss_res)
# ss_tot = np.sum((T2_avg_fit - np.mean(T2_avg_fit))**2)
# r_squared = 1 - (ss_res / ss_tot)
# print('R squared value is: ', r_squared)
#
# plt.plot(temp_fit, T2_avg_fit, marker = 'o')
# plt.plot(x_model, y_model, 'g--',
#          label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

plt.title('T2 as temperature is raised')
plt.xlabel('Temperature (mK)')
plt.ylabel("T2 (sec)")
plt.legend()
plt.show()