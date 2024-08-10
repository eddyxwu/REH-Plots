import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from scipy import special
from scipy.optimize import curve_fit

Q2_T2 = np.array([70.47, 63.14, 58.04, 60.73, 53.45, 43.64, 45.32, 35.77, 29.38, 20.93, 13.86, 13.61, 8.574, 10.12, 6.946, 5.864, 7.16]) * 1e-6
Q3_T2 = np.array([56.43, 50.33, 42.23, 36.91, 41.01, 24.03, 22.13, 14.93, 14.62, 8.404, 10.53, 7.709, 6.003, 4.743]) * 1e-6
Q4_T2 = np.array([48.96, 51.6, 46.14, 37.93, 32.72, 24.38, 20.98, 16.7, 14.79, 10.26, 10.57, 7.702, 7.041, 3.698]) * 1e-6

temp = np.array([13.79, 15.2, 36.45, 47.94, 50.77, 57.21, 69.97, 73.12, 88.12, 97.82, 114.85, 119.45, 138.25, 140.28, 156.03, 168.69, 184.98]) * 1e-3
temp_mod = np.array([15.2, 36.45, 50.77, 69.97, 73.12, 88.12, 97.82, 114.85, 119.45, 138.25, 140.28, 156.03, 168.69, 184.98]) * 1e-3

plt.plot(temp, Q2_T2, marker = 'o', color = 'b', label='Q2', alpha = 0.3)
# plt.plot(temp_mod, Q3_T2, marker = 'o', color='red', label='Q3', alpha = 0.3)
# plt.plot(temp_mod, Q4_T2, marker = 'o', color='orange', label='Q4', alpha = 0.3)

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
al_dk = 5.447400321 * (10 ** -23)
Q2_w = 3.0977 * 1e9 * 2 * pi
Q3_w = 3.175 * 1e9 * 2 * pi
Q4_w = 3.179 * 1e9 * 2 * pi

# def T2(temp, ne, dk, constant, kappa, w = Q2_w, chi = chi2):
#     x = (hbar * w) / (2 * k * temp)
#     Q = 1 / ((ne / pi) * (np.sqrt((2 * dk) / (hbar * w))) + (2 / pi) * np.exp(-dk / (k * temp)) * np.exp(x) * special.kn(0, x) * (1 + np.exp(-2 * x)))
#     T1 = Q / w
#
#     n = 1 / (np.exp((hbar * Q2_wc) / (k * temp) - 1))
#     gamma_phi = (kappa / 2) * (np.real(np.sqrt((1 + ((2 * 1j * chi) / kappa)) ** 2 + ((8 * 1j * (chi / kappa) * n))) - 1)) + constant
#     Tphi = 1 / gamma_phi
#
#     gamma_2 = (1 / Tphi) + (1 / (2 * T1))
#     return 1 / gamma_2
#
# popt, pcov = curve_fit(T2, temp, Q2_T2, p0 = [3.44 * 1e-7, 3 * 1e-23, 1.1 * 1e4, 3.68 * 1e5])
# ne_guess, dk_guess, constant_guess, kappa_guess = popt
# x_model = np.linspace(min(temp), max(temp), 100)
# y_model = T2(x_model, ne_guess, dk_guess, constant_guess, kappa_guess)
# plt.plot(x_model, y_model)
#
# plt.plot(x_model, y_model, 'b',
#          label='fit: ne=%5.3g, dk=%5.3g, constant=%5.3g, kappa=%5.3g' % tuple(popt))

# dk defined
# def T2(temp, ne, kappa, constant, dk = 3.2 * 1e-23, w = Q2_w, chi = chi2):
#     x = (hbar * w) / (2 * k * temp)
#     Q = 1 / ((ne / pi) * (np.sqrt((2 * dk) / (hbar * w))) + (2 / pi) * np.exp(-dk / (k * temp)) * np.exp(x) * special.kn(0, x) * (1 + np.exp(-2 * x)))
#     T1 = Q / w
#
#     n = 1 / (np.exp((hbar * Q2_wc) / (k * temp) - 1))
#     gamma_phi = (kappa / 2) * (np.real(np.sqrt((1 + ((2 * 1j * chi) / kappa)) ** 2 + ((8 * 1j * (chi / kappa) * n))) - 1)) + constant
#     Tphi = 1 / gamma_phi
#
#     gamma_2 = (1 / Tphi) + (1 / (2 * T1))
#     return 1 / gamma_2
#
# popt, pcov = curve_fit(T2, temp, Q2_T2, p0 = [0.5, kappa2, 9.96 * 1e3], maxfev=10000)
# ne_guess, kappa_guess, constant_guess = popt
# x_model = np.linspace(min(temp), max(temp), 100)
# y_model = T2(x_model, ne_guess, kappa_guess, constant_guess)
# plt.plot(x_model, y_model)
#
# plt.plot(x_model, y_model, 'b',
#          label='fit: ne=%5.3g, kappa=%5.3g, constant=%5.3g' % tuple(popt))

# remove constant
# def T2(temp, ne, kappa, dk = 3.2 * 1e-23, w = Q2_w, chi = chi2):
#     x = (hbar * w) / (2 * k * temp)
#     Q = 1 / ((ne / pi) * (np.sqrt((2 * dk) / (hbar * w))) + (2 / pi) * np.exp(-dk / (k * temp)) * np.exp(x) * special.kn(0, x) * (1 + np.exp(-2 * x)))
#     T1 = Q / w
#
#     n = 1 / (np.exp((hbar * Q2_wc) / (k * temp) - 1))
#     gamma_phi = (kappa / 2) * (np.real(np.sqrt((1 + ((2 * 1j * chi) / kappa)) ** 2 + ((8 * 1j * (chi / kappa) * n))) - 1))
#     Tphi = 1 / gamma_phi
#
#     gamma_2 = (1 / Tphi) + (1 / (2 * T1))
#     return 1 / gamma_2
#
# popt, pcov = curve_fit(T2, temp, Q2_T2, p0 = [100, kappa2], maxfev=1000000)
# ne_guess, kappa_guess = popt
# x_model = np.linspace(min(temp), max(temp), 100)
# y_model = T2(x_model, ne_guess, kappa_guess)
# plt.plot(x_model, y_model)
#
# plt.plot(x_model, y_model, 'b',
#          label='fit: ne=%5.3g, kappa=%5.3g' % tuple(popt))

# parameter check
def T2(temp, ne, dk, constant, kappa, w = Q2_w, chi = chi2):
    x = (hbar * w) / (2 * k * temp)
    Q = 1 / ((ne / pi) * (np.sqrt((2 * dk) / (hbar * w))) + (2 / pi) * np.exp(-dk / (k * temp)) * np.exp(x) * special.kn(0, x) * (1 + np.exp(-2 * x)))
    T1 = Q / w

    n = 1 / (np.exp((hbar * Q2_wc) / (k * temp) - 1))
    gamma_phi = (kappa / 2) * (np.real(np.sqrt((1 + ((2 * 1j * chi) / kappa)) ** 2 + ((8 * 1j * (chi / kappa) * n))) - 1)) + constant
    Tphi = 1 / gamma_phi

    gamma_2 = (1 / Tphi) + (1 / (2 * T1))
    return 1 / gamma_2

x_model = np.linspace(min(temp), max(temp), 100)
y_model = T2(x_model, 3.44 * 1e-7, 3 * 1e-23, 1.1 * 1e4, 3.68 * 1e5)
plt.plot(x_model, y_model, 'blue')

# def T2(temp, ne, dk, kappa, constant, constant2, w = Q2_w, chi = chi2):
#     x = (hbar * w) / (2 * k * temp)
#     Q = 1 / ((ne / pi) * (np.sqrt((2 * dk) / (hbar * w))) + (2 / pi) * np.exp(-dk / (k * temp)) * np.exp(x) * special.kn(0, x) * (1 + np.exp(-2 * x)))
#     T1 = Q / w
#
#     n = 1 / (np.exp((hbar * Q2_wc) / (k * temp) - 1))
#     gamma_phi = (kappa / 2) * (np.real(np.sqrt((1 + ((2 * 1j * chi) / kappa)) ** 2 + ((8 * 1j * (chi / kappa) * n))) - 1)) + constant
#     Tphi = 1 / gamma_phi
#
#     T2 = (1 / Tphi) + (1 / (2 * T1)) + constant2
#     return T2
#
# popt, pcov = curve_fit(T2, temp, Q2_T2, p0 = [0.5, al_dk, 0.5, 1, 100], maxfev=10000)
# ne_guess, dk_guess, constant_guess, constant2_guess, kappa_guess = popt
# x_model = np.linspace(min(temp), max(temp), 100)
# y_model = T2(x_model, ne_guess, dk_guess, constant_guess, constant2_guess, kappa_guess)
# plt.plot(x_model, y_model)
#
# plt.plot(x_model, y_model, 'b',
#          label='fit: ne=%5.3g, dk=%5.3g, constant=%5.3g, constant2=%5.3g, kappa=%5.3g' % tuple(popt))

plt.title('T2 as temperature is raised')
plt.xlabel('Temperature (mK)')
plt.ylabel("T2 (sec)")

plt.legend()
plt.show()


