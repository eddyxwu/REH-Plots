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

plt.plot(temp, Q2_Tphi, marker = 'o', color = 'b', label='Q2', alpha = 0.3)
plt.plot(temp_mod, Q3_Tphi, marker = 'o', color='red', label='Q3', alpha = 0.3)
plt.plot(temp_mod, Q4_Tphi, marker = 'o', color='orange', label='Q4', alpha = 0.3)

# model fit

# constants
pi = constants.pi
hbar = constants.hbar
k = constants.k
Q2_wc = 7.2026 * 1e9 * 2 * pi
Q3_wc = 7.2076 * 1e9 * 2 * pi
Q4_wc = 7.2286 * 1e9 * 2 * pi
chi2 = 252 * 1e3 * 2 * pi
chi3 = 263 * 1e3 * 2 * pi
chi4 = 264 * 1e3 * 2 * pi
kappa2 = 540 * 1e3 * 2 * pi
kappa3 = 220* 1e3 * 2 * pi
kappa4 = 320 * 1e3 * 2 * pi

def Generalization(temp, constant, kappa = kappa2, chi = chi2):
    n = 1 / (np.exp((hbar * Q2_wc) / (k * temp) - 1))
    gamma = (kappa / 2) * (np.real(np.sqrt((1 + ((2 * 1j * chi) / kappa))**2 + ((8 * 1j * (chi / kappa) * n))) - 1)) + constant
    return 1 / gamma

def Generalization_K(temp, constant, kappa):
    n = 1 / (np.exp((hbar * Q2_wc) / (k * temp) - 1))
    gamma = (kappa / 2) * (np.real(np.sqrt((1 + ((2 * 1j * chi2) / kappa))**2 + ((8 * 1j * (chi2 / kappa) * n))) - 1)) + constant
    return 1 / gamma

def Generalization_chi(temp, constant, chi):
    n = 1 / (np.exp((hbar * Q2_wc) / (k * temp) - 1))
    gamma = (kappa2 / 2) * (np.real(np.sqrt((1 + ((2 * 1j * chi) / kappa2))**2 + ((8 * 1j * (chi / kappa2) * n))) - 1)) + constant
    return 1 / gamma

def Generalization3(temp_mod, constant, kappa = kappa3, chi = chi3):
    n = 1 / (np.exp((hbar * Q2_wc) / (k * temp_mod) - 1))
    gamma = (kappa / 2) * (np.real(np.sqrt((1 + ((2 * 1j * chi) / kappa))**2 + ((8 * 1j * (chi / kappa) * n))) - 1)) + constant
    return 1 / gamma

def Generalization4(temp_mod, constant, kappa = kappa4, chi = chi4):
    n = 1 / (np.exp((hbar * Q2_wc) / (k * temp_mod) - 1))
    gamma = (kappa / 2) * (np.real(np.sqrt((1 + ((2 * 1j * chi) / kappa))**2 + ((8 * 1j * (chi / kappa) * n))) - 1)) + constant
    return 1 / gamma

def Photon(temp, constant, wc = Q2_wc):
    n = 1 / (np.exp((hbar * wc) / (k * temp) - 1))
    PureDephasing = n * (n + 1) * (chi2 ** 2 / kappa2) + constant
    return 1 / PureDephasing

def Photon3(temp, constant, chi, wc = Q3_wc):
    n = 1 / (np.exp((hbar * wc) / (k * temp) - 1))
    PureDephasing = n * (n + 1) * (chi ** 2 / kappa3) + constant
    return 1 / PureDephasing

def Photon4(temp, constant, chi, wc = Q4_wc):
    n = 1 / (np.exp((hbar * wc) / (k * temp) - 1))
    PureDephasing = n * (n + 1) * (chi ** 2 / kappa4) + constant
    return 1 / PureDephasing

# popt, pcov = curve_fit(Photon, temp, Q2_Tphi, p0 = [0.5])
# constant_guess = popt
# x_model = np.linspace(min(temp), max(temp), 100)
# y_model = Photon(x_model, constant_guess)
# plt.plot(x_model, y_model)
#
# plt.plot(x_model, y_model, 'b',
#          label='fit: constant=%5.3g' % tuple(popt))

# generalization
# x_model = np.linspace(min(temp), max(temp), 100)
# y_model = Generalization(x_model)
# plt.plot(x_model, y_model, 'g')

popt, pcov = curve_fit(Generalization, temp, Q2_Tphi, p0 = [0.5])
constant_guess = popt
x_model = np.linspace(min(temp), max(temp), 100)
y_model = Generalization(x_model, constant_guess)
plt.plot(x_model, y_model)

plt.plot(x_model, y_model, 'b',
         label='Q2 fit with set params: constant=%5.3g' % tuple(popt))

popt, pcov = curve_fit(Generalization3, temp_mod, Q3_Tphi, p0 = [0.5])
constant_guess = popt
x_model = np.linspace(min(temp_mod), max(temp_mod), 100)
y_model = Generalization3(x_model, constant_guess)
plt.plot(x_model, y_model)

plt.plot(x_model, y_model, 'r',
         label='Q3 fit with set params: constant=%5.3g' % tuple(popt))

popt, pcov = curve_fit(Generalization4, temp_mod, Q4_Tphi, p0 = [0.5])
constant_guess = popt
x_model = np.linspace(min(temp_mod), max(temp_mod), 100)
y_model = Generalization4(x_model, constant_guess)
plt.plot(x_model, y_model)

plt.plot(x_model, y_model, 'orange',
         label='Q4 fit with set params: constant=%5.3g' % tuple(popt))

popt, pcov = curve_fit(Generalization_K, temp, Q2_Tphi, p0 = [0.5, 100])
constant_guess, kappa_guess = popt
x_model = np.linspace(min(temp), max(temp), 100)
y_model = Generalization_K(x_model, constant_guess, kappa_guess)
plt.plot(x_model, y_model)

plt.plot(x_model, y_model, 'black',
         label='fit: constant=%5.3g, kappa=%5.3g' % tuple(popt))
#
# popt, pcov = curve_fit(Generalization_chi, temp, Q2_Tphi, p0 = [0.5, 100])
# constant_guess, chi_guess = popt
# x_model = np.linspace(min(temp), max(temp), 100)
# y_model = Generalization_chi(x_model, constant_guess, chi_guess)
# plt.plot(x_model, y_model)
#
# plt.plot(x_model, y_model, 'green',
#          label='fit: constant=%5.3g, chi=%5.3g' % tuple(popt))

# popt, pcov = curve_fit(Photon3, temp_mod, Q3_Tphi, p0 = [1, 400])
# constant_guess, chi_guess = popt
# x_model = np.linspace(min(temp), max(temp), 100)
# y_model = Photon3(x_model, constant_guess, chi_guess)
# plt.plot(x_model, y_model)
#
# plt.plot(x_model, y_model, 'r',
#          label='fit: constant=%5.3g, chi=%5.3g' % tuple(popt))
#
# popt, pcov = curve_fit(Photon4, temp_mod, Q4_Tphi, p0 = [1, 400])
# constant_guess, chi_guess = popt
# x_model = np.linspace(min(temp), max(temp), 100)
# y_model = Photon4(x_model, constant_guess, chi_guess)
# plt.plot(x_model, y_model)
#
# plt.plot(x_model, y_model, 'orange',
#          label='fit: constant=%5.3g, chi=%5.3g' % tuple(popt))

plt.title('T_phi as temperature is raised')
plt.xlabel('Temperature (K)')
plt.ylabel("Time (sec)")

plt.legend()
plt.show()

# test_times = Generalization(temp, 9.96 * 1e3, kappa2, chi2)
# print(test_times)