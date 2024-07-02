import numpy as np
from os.path import dirname, abspath
import matplotlib.pyplot as plt

def best_fit_line_plot(x, y, x_label, y_label, title):
    '''
    Plot data in scatterplot with best fit line
    '''
    slope, intercept = np.polyfit(x, y, 1)
    best_fit_line = slope * x+ intercept
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(x, y, color='blue', label='Data points')
    ax1.plot(x, best_fit_line, color='red', label='Best fit line')
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    ax1.set_title(title)


d = dirname(dirname(dirname(dirname(abspath(__file__)))))
if __name__ == "__main__":
    filepath = d+'/Europa/EuropaProfile_CustomSolution_0.0ppt_Tb268.305K.txt'
    data = np.genfromtxt(filepath, skip_header = 82)
    data_under_140_MPa = data[(data[:, 0] < 140) & (data[:, 0] > 50)]
    pressures_under_140_MPa = data_under_140_MPa[:, 0].T
    densities = data_under_350_MPa = data_under_140_MPa[:, 4].T
    best_fit_line_plot(pressures_under_140_MPa, densities, "Pressure (MPa)", "Density (kg/m3)", "Density vs Pressure")
    plt.show()
