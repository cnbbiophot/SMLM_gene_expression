#Funciones genéricas que usa anaCellC
#Pueden ser usadas en otros programas
#Funciones:

# Funciones de ajuste de anaCellC: línea recta, exponencial, acumulativas, etc
#
# jri - 31.1.24


from os.path import isfile, join
import numpy as np
from scipy import optimize, stats
import re #Para modificar el archivo param.dat


def fit_cumFrame_linear (xData, yData):
    #Ajusta la acumulativa de los knock out. 
    #x0 es el guess en formato m, y0 (pendiente y ordenada en el origen)
    x0=(0.2, 1)
    yErr=yData**0.5
    yErr[np.where(yErr==0)]=1
    paramFit = optimize.least_squares(err_straightline, x0, args=(xData, yData, yErr), method='lm')
    return paramFit

def fit_cumFrame_exp (xData, yData, yErr, x0=None):
    #Ajusta la acumulativa a 1 - exp
    # y = A*(1-exp(-frame/tau))
    #x0 es el guess (A, tau). 
    if x0 is None:
        x0=(100, 5)
    yErr[np.where(yErr==0)]=1
    xData=xData.astype(np.float64)
    yData=yData.astype(np.float64)
    paramFit = optimize.least_squares(err_exp_complementaria, x0, bounds=((0, 0), (np.inf, np.inf)), args=(xData, yData, yErr))
    return paramFit

def fit_cumFrame_gamma (xData, yData, yErr, x0=None):
    #Ajusta la acumulativa a una cdf de gamma. 
    # Sería Erlang si hubiese sido una resta de dos exponenciales con el mismo tau
    #Pero como son distintas, entonces es gamma
    #x0 es el guess (A, forma, tau). 
    if x0 is None:
        x0=(2.2, 2, 20)
    yErr[np.where(yErr==0)]=1
    xData=xData.astype(np.float64)
    yData=yData.astype(np.float64)
    paramFit = optimize.least_squares(err_gamma_cdf, x0, bounds=((0, 0, 0), (np.inf, np.inf, np.inf)), args=(xData, yData, yErr))
    return paramFit


def fit_cumFrame_exp_linear (xData, yData, yErr=None):
    #Ajusta la acumulativa de de las no knock out a un logaritmo y una recta
    # y = A*(1-exp(-frame/tau))+m*frame
    #x0 es el guess (A, tau, m). 
    x0=(6, 5, 0.003)
    if yErr is None:
        yErr=yData**0.5
    yErr[np.where(yErr==0)]=1
    xData=xData.astype(np.float64)
    yData=yData.astype(np.float64)
    paramFit = optimize.least_squares(err_exp_lineal_complementaria, x0, bounds=((0, 0, 0), (np.inf, np.inf, np.inf)), args=(xData, yData, yErr))
    return paramFit

def fit_cumFrame_exp_linear_2 (tau, xData, yData, yErr=None):
    #Ajusta la acumulativa de de las no knock out a un logaritmo y una recta
    # y = A*(1-exp(-frame/tau))+m*frame
    #x0 es el guess (A, m). 
    x0=(6, 0.003)
    if yErr is None:
        yErr=yData**0.5
    yErr[np.where(yErr==0)]=1
    xData=xData.astype(np.float64)
    yData=yData.astype(np.float64)
    paramFit = optimize.least_squares(
        err_exp_lineal_complementaria_2, x0, bounds=((0, 0), (np.inf, np.inf)), args=(tau, xData, yData, yErr))
    return paramFit


def fit_doble_exp_complementaria (A_KO, tau_KO, xData, yData, yErr=None):
    #Ajusta la acumulativa de de las no knock out a un logaritmo y una recta
    # y = A*(1-exp(-frame/tau))+A_KO*(1-exp(-frame/tau_KO))+m*frame
    # A_KO y tau_KO fijados
    # x0 es el guess (A, tau, m). 
    x0=(6, 50, 0.039)
    if yErr is None:
        yErr=yData**0.5
    yErr[np.where(yErr==0)]=1
    xData=xData.astype(np.float64)
    yData=yData.astype(np.float64)
    # paramFit = \
    # optimize.least_squares(err_doble_exp_complementaria, x0, bounds=((0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf)), args=(tau_KO, xData, yData, yErr))
    paramFit = \
    optimize.least_squares(err_doble_exp_complementaria, x0, bounds=((0, 0, 0), (np.inf, np.inf, np.inf)), args=(A_KO, tau_KO, xData, yData, yErr))

    return paramFit

def fit_doble_exp_complementaria_libre (xData, yData, yErr=None):
    #Ajusta la acumulativa de de las no knock out a la cdf exponencial y una recta
    # y = A*(1-exp(-frame/tau))+A_KO*(1-exp(-frame/tau_KO))+m*frame
    # TODOS LOS PARÁMETROS LIBRES
    # x0 es el guess (A, tau, A_KO, tau, m). 
    x0=(6, 7, 3, 50, 0.039)
    if yErr is None:
        yErr=yData**0.5
    yErr[np.where(yErr==0)]=1
    xData=xData.astype(np.float64)
    yData=yData.astype(np.float64)
    paramFit = optimize.least_squares(err_doble_exp_complementaria_libre, x0,
        bounds=((0, 0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf, np.inf)), args=(xData, yData, yErr))

    return paramFit


def fit_expDecay (xData, yData):
    # y = A*exp(-frame/tau)+B
    #x0 es el guess (A, tau, B). 
    x0=(yData[0], 5, yData[-1])
    yErr=yData**0.5
    yErr[np.where(yErr==0)]=1
    xData=xData.astype(np.float64)
    yData=yData.astype(np.float64)
    paramFit = optimize.least_squares(err_expDecay, x0, bounds=((0, 0, 0), (np.inf, np.inf, np.inf)), args=(xData, yData, yErr))
    return paramFit

def err_exp_lineal_complementaria(x, xData, yData, yErr):
    yFit=exp_lineal_complementaria(x, xData)
    err=(yData-yFit)/yErr
    return err

def err_exp_complementaria(x, xData, yData, yErr):
    yFit=exp_complementaria(x, xData)
    err=(yData-yFit)/yErr
    return err

def exp_complementaria (params, xData):
    A=params[0]
    tau=params[1]
    y=A*(1-np.exp(-xData/tau))
    return y

def exp_lineal_complementaria (params, xData):
    A=params[0]
    tau=params[1]
    m=params[2]
    y=A*(1-np.exp(-xData/tau))+m*xData
    return y

def err_exp_lineal_complementaria_2 (x, tau, xData, yData, yErr):
    # tau es un parámetro fijado
    yFit=exp_lineal_complementaria((x[0], tau, x[1]), xData)
    err=(yData-yFit)/yErr
    return err

def err_doble_exp_complementaria (x, A_KO, tau_KO, xData, yData, yErr):
    # A_KOy tau_KO son  parámetros fijados
    # yFit=doble_exp_complementaria((x[0], x[1], x[2], x[3]), tau_KO, xData)
    yFit=doble_exp_complementaria((x[0], x[1], A_KO, tau_KO, x[2]), xData)
    err=(yData-yFit)/yErr
    return err

def err_doble_exp_complementaria_libre (x, xData, yData, yErr):
    yFit=doble_exp_complementaria(x, xData)
    err=(yData-yFit)/yErr
    return err


def err_straightline (x, xData, yData, yErr):
    m=x[0]
    y0=x[1]
    #y0=yData[0] #Fijo el y0 al primer dato de yData
    yFit=straightline((m, y0), xData)
    err=(yData-yFit)/yErr
    return err

def err_expDecay (x, xData, yData, yErr):
    yFit=expDecay(x, xData)
    err=(yData-yFit)/yErr
    return err

def err_gamma_cdf (x, xData, yData, yErr):
    yFit=gamma_cdf(x, xData)
    err=(yData-yFit)/yErr
    return err

def straightline(params, xData):
    m=params[0]
    y0=params[1]
    y=m*xData+y0
    return y

def expDecay(params, xData):
    I0=params[0]
    tau=params[1]
    B=params[2]
    y=I0*np.exp(-(xData-xData[0])/tau)+B
    return y

def gamma_cdf(params, xData):
    A=params[0]
    shape=params[1] #Factor de forma
    tau=params[2] # factor de escala
    y=A*stats.gamma.cdf(xData, shape, 0, tau)
    return y


def doble_exp_complementaria (params, xData):
    A=params[0]
    tau=params[1]
    A_KO=params[2]
    tau_KO=params[3]
    m=params[4]
    y=A*(1-np.exp(-xData/tau))+A_KO*(1-np.exp(-xData/tau_KO))+m*xData
    return y
