import matplotlib.pyplot as plt
import numpy as np

power = np.array([0, 0, 30, 50, 60, 75, 115, 125, 180, 225, 300, 337.5, 420, 450, 520, 635, 709])
temp = np.array([13.79, 15.2, 36.45,47.94,50.77,57.21,69.97,73.12,88.12,97.82,114.85,119.09,138.25, 140.28, 156.03, 168.69, 184.98])

plt.plot(power, temp, marker = 'o', color = 'purple')
plt.title('Temperature Dependence with Power')
plt.ylabel('Temperature (mK)')
plt.xlabel('Power (μW)')

plt.show()
