import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from scipy import special
from scipy.optimize import curve_fit



Q2_T1 = np.array([119.4, 107.8, 81.58, 95.81, 108.9, 53.69, 94, 64.8, 78.13, 75.18, 79.83, 81.28, 52.45, 68.94, 27.28, 25.25, 10.32]) * 1e-6
Q3_T1 = np.array([91.23, 63.86, 66.17, 89.79, 140.9, 59.3, 93.57, 79.84, 61.12, 73.23, 80.54, 44.51, 24.53, 6.927]) * 1e-6
Q4_T1 = np.array([69.58, 109.5, 97.27, 92.64, 89.97, 80.01, 61.74, 63.74, 107.4, 55.03, 64.86, 36.8, 24.03, 12]) * 1e-6

temp = np.array([13.79, 15.2, 36.45, 47.94, 50.77, 57.21, 69.97, 73.12, 88.12, 97.82, 114.85, 119.45, 138.25, 140.28, 156.03, 168.69, 184.98]) * 1e-3
temp_mod = np.array([15.2, 36.45, 50.77, 69.97, 73.12, 88.12, 97.82, 114.85, 119.45, 138.25, 140.28, 156.03, 168.69, 184.98]) * 1e-3

plt.plot(temp, Q2_T1, marker = 'o', color = 'b', label='Q2', alpha = 0.3)
plt.plot(temp_mod, Q3_T1, marker = 'o', color='red', label='Q3', alpha = 0.3)
plt.plot(temp_mod, Q4_T1, marker = 'o', color='orange', label='Q4', alpha = 0.3)

# model fit

# constants
pi = constants.pi
hbar = constants.hbar
k = constants.k
al_dk = 5.447400321 * (10 ** -23)
Q2_w = 3.0977 * 1e9 * 2 * pi
Q3_w = 3.175 * 1e9 * 2 * pi
Q4_w = 3.179 * 1e9 * 2 * pi

def QP(temp, ne, dk, w = Q2_w):
    x = (hbar * w) / (2 * k * temp)
    Q = 1 / ((ne / pi) * (np.sqrt((2 * dk) / (hbar * w))) + (2 / pi) * np.exp(-dk / (k * temp)) * np.exp(x) * special.kn(0, x) * (1 + np.exp(-2 * x)))
    return Q / w

def QP_Q3(temp, ne, dk, w = Q3_w):
    x = (hbar * w) / (2 * k * temp)
    Q = 1 / ((ne / pi) * (np.sqrt((2 * dk) / (hbar * w))) + (2 / pi) * np.exp(-dk / (k * temp)) * np.exp(x) * special.kn(0, x) * (1 + np.exp(-2 * x)))
    return Q / w

def QP_Q4(temp, ne, dk, w = Q4_w):
    x = (hbar * w) / (2 * k * temp)
    Q = 1 / ((ne / pi) * (np.sqrt((2 * dk) / (hbar * w))) + (2 / pi) * np.exp(-dk / (k * temp)) * np.exp(x) * special.kn(0, x) * (1 + np.exp(-2 * x)))
    return Q / w

popt, pcov = curve_fit(QP, temp, Q2_T1, p0 = [0.5, al_dk])
popt
ne, dk_guess = popt
x_model = np.linspace(min(temp), max(temp)*1.5, 100)
y_model = QP(x_model, ne, dk_guess)
plt.plot(x_model, y_model)

plt.plot(x_model, y_model, 'b',
         label='fit: ne=%5.3g, dk=%5.3g' % tuple(popt))

popt, pcov = curve_fit(QP_Q3, temp_mod, Q3_T1, p0 = [0.5, al_dk])
popt
ne, dk_guess = popt
x_model = np.linspace(min(temp), max(temp)*1.5, 100)
y_model = QP_Q3(x_model, ne, dk_guess)
plt.plot(x_model, y_model)

plt.plot(x_model, y_model, 'r',
         label='fit: ne=%5.3g, dk=%5.3g' % tuple(popt))

popt, pcov = curve_fit(QP_Q4, temp_mod, Q4_T1, p0 = [0.5, al_dk])
popt
ne, dk_guess = popt
x_model = np.linspace(min(temp), max(temp)*1.5, 100)
y_model = QP_Q4(x_model, ne, dk_guess)
plt.plot(x_model, y_model)

plt.plot(x_model, y_model, 'orange',
         label='fit: ne=%5.3g, dk=%5.3g' % tuple(popt))

plt.title('T1 as temperature is raised')
plt.xlabel('Temperature (K)')
plt.ylabel("T1 (sec)")

plt.legend()
plt.show()

# test_times = QP(temp, 3.44 * 1e-7, 3.2 * 1e-23, w = Q2_w)
# print(test_times)