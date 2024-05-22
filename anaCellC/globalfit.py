#Global fitting
#Adaptado de : https://stackoverflow.com/questions/51482486/python-global-fitting-for-data-sets-of-different-sizes
# y de : https://stackoverflow.com/questions/59505194/global-fitting-using-scipy-curve-fit
#
# jri - 30Abr24

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import anaCellC.utils as u

# x1 = np.array([5.0, 6.1, 7.2, 8.3])
# y1 = np.array([ 16.00,  18.42,  20.84,  23.26])
# x2 = np.array([15.0, 16.1, 17.2, 18.3, 19.4])
# y2 = np.array([-20.00, -25.50, -31.00, -36.50, -42.00])

# # some initial parameter values
# initialParameters = np.array([1.0, 1.0, 1.0])

# def function1(data, a, b, c): # not all parameters are used here, c is shared
#         return a * data + c

# def function2(data, a, b, c): # not all parameters are used here, c is shared
#         return b * data + c

def function1(xData, A, tau, A_KO, tau_KO, m): # not all parameters are used here, A_KO, tau_KO and m are shared
        y = u.doble_exp_complementaria ((A, tau, A_KO, tau_KO, m), xData)
        return y

def function2(xData, A_KO, tau_KO, m): # not all parameters are used here, A_KO, tau_KO and m are shared
        y = u.exp_lineal_complementaria ((A_KO, tau_KO, m), xData)
        return y




def combinedfunction(comboX, A, tau, A_KO, tau_KO, m):
    # single data reference passed in, extract separate data
    #Esto lo hago así porque en mi caso sé que x1 y x2 tienen el mismo tamaño. 
    #En realidad tendría que pasar el tamaño como parámetro
    x1=comboX[:int(len(comboX)/2)]
    extract1 = comboX[:len(x1)] # first data
    extract2 = comboX[len(x1):] # second data

    result1 = function1(extract1, A, tau, A_KO, tau_KO, m)
    result2 = function2(extract2, A_KO, tau_KO, m)

    return np.append(result1, result2)


def globalfit (x1, y1,  yErr1, x2, y2, yErr2, initialParameters):
    comboY = np.append(y1, y2)
    comboX = np.append(x1, x2)
    comboErr = np.append(yErr1, yErr2)

    if len(y1) != len(x1):
        raise(Exception('Unequal x1 and y1 data length'))
    if len(y2) != len(x2):
        raise(Exception('Unequal x2 and y2 data length'))

    # curve fit the combined data to the combined function
    fittedParameters, pcov = curve_fit(combinedfunction, comboX, comboY, p0=initialParameters, sigma=comboErr, absolute_sigma=True)
    # fittedParameters, pcov = curve_fit(combinedFunction, comboX, comboY, p0=initialParameters)

    # values for display of fitted function
    A, tau, A_KO, tau_KO, m = fittedParameters

    y_fit_1 = function1(x1, A, tau, A_KO, tau_KO, m) # first data set, first equation
    y_fit_2 = function2(x2, A_KO, tau_KO, m) # second data set, second equation

    # fig,ax = plt.subplots(1, 1)
    # ax.plot(x1, y1) # plot the raw data
    # ax.plot(x2, y2) # plot the raw data
    # ax.plot(x1, y_fit_1) # plot the equation using the fitted parameters
    # ax.plot(x2, y_fit_2) # plot the equation using the fitted parameters
    # # plt.show()

    # print('A, tau, A_KO, tau_KO, m: ', fittedParameters)
    return A, tau, A_KO, tau_KO, m, y_fit_1, y_fit_2

# globalfit (x1, y1,  x2, y2, initialParameters)