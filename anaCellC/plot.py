# Funciones de plot de anaCellC
#
# Las importantes on histo_x1 y histocumulative_x1 y calculaerrorhisto
#
# plot2Graphs: #Representa el recuento acumulado de moléculas de células y el total
# idx2frame
# plotHisto
# plotCountsPerFrame
# plot3Graphs
# histocumulative_x2: Representa dos histogramas superpuestos junto con sus distribuciones acumuladas
# histocumulative_x1: Representa un histograma (sólo uno) junto con su distribución acumulada
# cambiacolorbarrashisto: Cambia el color de las barras de un histograma
# histo_x1: Representa un histograma junto con sus errores. Sin CDF
# calculaerrorhisto: Calcula los errores en los histogramas de frecuencia
#
# jri - 23Feb24

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def plot2Graphs (numContornos, cumCell, cumFrame, yFit, frameFit, frameQuant):
    #Representa las localizaciones por célula y acumuladas en el frame
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('axes', labelsize=14)
    plt.rc('legend', fontsize=12)

    # Representa el acumulado en el frame
    fig_cum, (ax_cell, ax_total) = plt.subplots(2, 1, sharex=True)
    fig_cum.subplots_adjust(hspace=0)
    ax_total.plot(cumFrame[:, 0], cumFrame[:, 1], label='Count inside cells')

    #Representa las localizaciones por célula. El área control es siempre la última (con índice numContornos-1)
    k = 0
    while k < numContornos - 1:  # Todas menos el área de control. k es un índice que empieza en 0 y acaba en (numContornos-2) incluido. 
    #Por tanto cuenta todas las cúlulas pero no el área control
        idx = (cumCell[:, 0] == k)  # Encuentra la célula correspondiente
        #Convierto los índices discretos en valores de todos los frames:
        cumCell_cont=idx2frame(cumCell[idx, 2], cumCell[idx, 3], frameFit)
        ax_cell.plot(cumCell_cont[:, 0], cumCell_cont[:, 1], label=str(k))
        #Representación anterior:
        #ax_cell.plot(cumCell[idx, 2], cumCell[idx, 3], label=str(k))
        k = k + 1

    #Representa el ajuste y el intervalo de cuantificación
    if len(yFit)>0:
        ax_total.plot(cumFrame[:, 0], yFit, label='Model')
        #Representa el límite del intervalo de cuantificación
        ylim=ax_total.get_ylim()
        ax_total.plot([frameQuant, frameQuant], ylim, linestyle='dashed', linewidth=1)
    '''
    #Añado para que representa también el área control multiplicada por el número de contornos
    idx = (cumCell[:, 0] == numContornos-1)  # Encuentra el área de control
    cumCell_cont=idx2frame(cumCell[idx, 2], cumCell[idx, 3], frameFit)
    ax_total.plot(cumCell_cont[:,0], cumCell_cont[:,1]*(numContornos-1), label='Control area')
    '''

    ax_cell.set_ylabel('Fluor. events')
    ax_total.set_xlabel('Frame')
    ax_total.set_ylabel('Fluor. events')
    #ax_cell.set_title("Cumulative molecule count")
    #ax_cell.set_xscale('log')
    leg_cell = ax_cell.legend(loc='upper right')
    leg_total = ax_total.legend()
    plt.pause(0.03)
    return fig_cum

def plotcumframe (numContornos, cumFrame, yFit, frameFit, frameQuant):
    #Representa las localizaciones acumuladas en el frame
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('axes', labelsize=14)
    plt.rc('legend', fontsize=12)

    # Representa el acumulado en el frame
    fig_cum, ax_total = plt.subplots(1, 1)
    ax_total.plot(cumFrame[:, 0], cumFrame[:, 1], label='Fluor. events')

    #Representa el ajuste y el intervalo de cuantificación
    if len(yFit)>0:
        ax_total.plot(cumFrame[:, 0], yFit, label='Model')
        #Representa el límite del intervalo de cuantificación
        ylim=ax_total.get_ylim()
        ax_total.plot([frameQuant, frameQuant], ylim, linestyle='dashed', linewidth=1)

    ax_total.set_xlabel('Frame')
    ax_total.set_ylabel('Fluor. events')
    #ax_total.set_title("Cumulative molecule count")
    #ax_total.set_xscale('log')
    leg_total = ax_total.legend()
    plt.pause(0.03)
    return fig_cum




def idx2frame (fNum, count, frameFit):
    #fNum es el número del frame en el que ha habido (al menos) una localización. Count es el número de localizaciones en ese frame
    #fNum puede empezar en 1 o no
    #Lo convierto en cumCell_c
    frameFit=int(frameFit)
    cumCell_c=np.zeros((frameFit, 2))
    cumCell_c[:,0]=np.arange(1,frameFit+1)
    fMax=len(fNum)
    if fMax>0: #Las células en las que no hay cuentas van vacías. 
        fNum=fNum.astype('int')
        ndx=0
        if fNum[0]>1:
            n0=0; nFin=fNum[0]-1 #Esto son índices
            cumCell_c[n0:nFin,1]=0
        while ndx<fMax-1: #Este bucle va de 0 a fMax-2, incluidos
            n0=fNum[ndx]-1; nFin=fNum[ndx+1]-1
            cumCell_c[n0:nFin,1]=count[ndx]
            ndx=ndx+1
        #Ahora hago el último índice:
        #Si fNum[fMax-1]==frameFit no es problema; si fNum[fMax-1]<frameFit tampoco es problema porque completa hasta frameMit
        n0=fNum[fMax-1]-1; nFin=frameFit
        cumCell_c[n0:nFin,1]=count[fMax-1]
        #numpy slice: Me quedo con las fMax primeras filas
        cumCell_c=cumCell_c[:fNum[fMax-1],:]
    else:
        cumCell_c=np.zeros((2, 2))
    return cumCell_c


def plotHisto (numLocaCell):
    fig_hist, ax_hist = plt.subplots(1, 1)
    # Quito el área control (la última célula)
    counts, bins = np.histogram(numLocaCell[:-1, 1], bins=15)
    ax_hist.stairs(counts, bins, hatch='//')
    ax_hist.set_title("Histogram of fluor. events per cell")
    ax_hist.set_xlabel('Fluorescent events')
    ax_hist.set_ylabel('Cell count')
    plt.pause(0.03)

    return fig_hist

    #plt.show(block=False)
    #plt.show()
    #show() bloquea la ejecución.

def plotCountsPerFrame (countsPerFrame, yFit=[], countsControlArea=None):
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('axes', labelsize=14)
    plt.rc('legend', fontsize=12)
    fig_frame, ax_total = plt.subplots(1, 1)
    
    numContornos=np.shape(countsPerFrame)[1]-1 #La primera columna es el índice del frame
    frame=countsPerFrame[:,0]
    if numContornos==1:
        ax_total.plot(frame, countsPerFrame[:,1], label='Count inside cells')
    else:
        ax_total.plot(frame, np.sum(countsPerFrame[:,1:-1], axis=1), label='Count inside cells')

    #Representa el ajuste y el intervalo de cuantificación
    if len(yFit)>0:
        ax_total.plot(frame, yFit, label='Fit')
    #     #Representa el límite del intervalo de cuantificación
    #     ylim=ax_total.get_ylim()
    #     ax_total.plot([frameQuant, frameQuant], ylim, linestyle='dashed', linewidth=1)
    if countsControlArea is not None:
        ax_total.plot(countsControlArea[:,0], countsControlArea[:,1], label='Control')

    ax_total.set_xlabel('Frame')
    ax_total.set_ylabel('Fluor. events')
    leg_total = ax_total.legend()
    plt.pause(0.03)
    return fig_frame, ax_total   
    

def plot3Graphs (numContornos, numLocaCell, cumCell, cumFrame, paramFit, filePath, rootName):
    # fig=plt.figure();
    #ax_cell= fig.add_subplot(311)
    # ax_noise=fig.add_subplot(312)
    #ax_total = fig.add_subplot(313)

    fig_cum, (ax_cell, ax_noise, ax_total) = plt.subplots(3, 1, sharex=True)
    fig_cum.subplots_adjust(hspace=0)
    # Encuentra el área de control
    idx = (cumCell[:, 0] == numContornos - 1)
    ax_noise.plot(cumCell[idx, 2], cumCell[idx, 3], label='Control area')

    ax_total.plot(cumFrame[:, 0], cumFrame[:, 1], label='Total mol. inside cells')
    #yFit=paramFit.x[0]*cumFrame[:, 0]+paramFit.x[1]
    yFit=u.exp_complementaria (paramFit.x, cumFrame[:, 0].astype(np.float64))
    ax_total.plot(cumFrame[:, 0], yFit)
    k = 0
    while k < numContornos - 2:  # Todas menos el área de control
        idx = (cumCell[:, 0] == k)  # Encuentra la célula correspondiente
        ax_cell.plot(cumCell[idx, 2], cumCell[idx, 3], label=str(k))
        k = k + 1

    ax_cell.set_ylabel('Mol. count')
    ax_total.set_xlabel('Frame')
    ax_noise.set_ylabel('Mol. count')
    ax_total.set_ylabel('Mol. count')
    ax_cell.set_title("Cumulative molecule count")
    #ax_cell.set_xscale('log')
    leg_noise = ax_noise.legend()
    leg_cell = ax_cell.legend()
    leg_total = ax_total.legend()

    fig_hist, ax_hist = plt.subplots(1, 1)
    # Quito el área control (la última célula)
    counts, bins = np.histogram(numLocaCell[:-1, 1], bins=15)
    ax_hist.stairs(counts, bins, hatch='//')
    ax_hist.set_title("Molecule count histogram")
    ax_hist.set_xlabel('Number of molecules')
    ax_hist.set_ylabel('Cell count')

#        if saveFigs:
    print('Saving ' + rootName + '_cum.png')
    fig_cum.savefig(join(filePath, rootName + '_cum.png'))
    print('Saving' + rootName + '_hist.png')
    fig_hist.savefig(join(filePath, rootName + '_hist.png'))
    #plt.show(block=False)
    plt.show()
    #plt.pause(0.03)
    #show() bloquea la ejecución.
    
    #input ('Press any key...')
    #plt.close('all')

def histocumulative_x2 (data1, data2, num_bins=None, min_histo=0, max_histo=None, label1='Data 1', label2='Data 2'):
    #Representa dos histogramas superpuestos junto con su distribución acumulada
    #Lo uso para representar, por ejemplo, los resultados de KO y Ocre
    # data1 va en azul por debajo
    # data2 va en amarillo por encima
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('axes', labelsize=14)
    plt.rc('legend', fontsize=12)


    n1=len(data1)
    n2=len(data2)
    if num_bins == None:
        num_bins=np.sqrt(np.max([n1, n2]))
    if max_histo != None:
        bin_edges= np.linspace(min_histo, max_histo, num_bins+1)
        

    fig, (ax_hist, ax_cum) = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    #Tengo que usar esto para representar correctamente la densidad de probabilidad de en el caso de las pmf
    # freq1, _, bars1=ax_hist.hist(data1, bins=bins_histo, weights=np.ones(n1)/n1, density=False, rwidth=0.85, alpha=0.5, label=label1)
    # freq2, _, bars2=ax_hist.hist(data2, bins=bins_histo, weights=np.ones(n2)/n2,density=False, rwidth=0.85, alpha=0.5, label=label2)
    #Esto queda feo porque pone una línea al final
    # c_KO=ax_cum.hist(intensidad_KO, bins=np.linspace(0,lim_intensidad, n_KO+1),density='true', cumulative='True', histtype="step", alpha=0.5, label='KO')
    # c_ocre=ax_cum.hist(intensidad_ocre, bins=np.linspace(0,lim_intensidad, n_ocre+1), density='true', cumulative='True', histtype="step", alpha=0.5, label='ocre')
    #Para versiones nuevas de matplotlib uso
    #c_KO=ax_cum.ecdf(intensidad_KO)
    #c_ocre=ax_cum.ecdf(intensidad_ocre)
    
    '''
    #Usaría esto si quisiese hacer una suma acumulativo con bins de 1 en 1
    freq_sin_bin_1=np.histogram(data1, bins=np.linspace(0,max_histo, n1+1), density=False, weights=np.ones(n1)/n1)
    freq_sin_bin_2=np.histogram(data2, bins=np.linspace(0,max_histo, n2+1), density=False, weights=np.ones(n2)/n2)
    cum1=np.cumsum(freq_sin_bin_1[0])
    cum2=np.cumsum(freq_sin_bin_2[0])
    '''
    weights1=np.ones(n1)/n1
    weights2=np.ones(n2)/n2
    freq1, _=np.histogram(data1, bins=bin_edges, weights=weights1, density=False)
    counts1, *_ = np.histogram(data1, bins=bin_edges, density=False)
    freq2, _=np.histogram(data2, bins=bin_edges, weights=weights2, density=False)
    counts2, *_ = np.histogram(data2, bins=bin_edges, density=False)
    errors1, bin_centers=calculaerrorhisto(data1, weights1, bin_edges)
    errors2, bin_centers=calculaerrorhisto(data2, weights2, bin_edges)
    width=(bin_edges[0]-bin_edges[1])*0.85
    bars1= ax_hist.bar(bin_centers, freq1, width=width, yerr=errors1, label=label1, alpha=0.5)
    bars2= ax_hist.bar(bin_centers, freq2, width=width, yerr=errors2, label=label2, alpha=0.5)

    cum1=np.cumsum(freq1)
    cum2=np.cumsum(freq2)
    #Añado un punto al comienzo para hacer el primer step completo en el acumulativo y que empiece desde 0
    cum1=np.insert(cum1,0,0) 
    cum2=np.insert(cum2,0,0)
    plot_cum1=ax_cum.step(bin_edges, cum1, alpha=0.5, label=label1)
    plot_cum2=ax_cum.step(bin_edges, cum2, alpha=0.5, label=label2)
    leg=ax_hist.legend()
    #ax_hist.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.02e'))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(-2,3), useMathText=True)

    ylabel=ax_hist.set_ylabel('Frequency', fontsize=14)
    ylabel=ax_cum.set_ylabel('Cumulative frequency', fontsize=14)

    return fig, ax_hist, ax_cum, bars1, bars2, plot_cum1, plot_cum2, bin_centers, freq1, errors1, freq2, errors2, bin_edges, cum1, cum2

def histocumulative_x1 (data, num_bins=None, min_histo=0, max_histo=None, label='Data'):
    #Representa un histograma junto con su distribución acumulada
    #Lo uso para representar, por ejemplo, los resultados de ocre, ambar o WT por separado
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('axes', labelsize=14)
    plt.rc('legend', fontsize=12)

    n=len(data)
    if num_bins == None:
        num_bins=np.sqrt(np.max(n))
    if max_histo != None:
        bin_edges= np.linspace(min_histo, max_histo, num_bins+1)

    fig, (ax_hist, ax_cum) = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    #Tengo que usar esto para representar correctamente la densidad de probabilidad de en el caso de las pmf
    # freq, bins, bars=ax_hist.hist(data, bins=bins_histo, weights=np.ones(n)/n, density=False, rwidth=0.85, alpha=1, label=label)

    freq, *_= np.histogram(data, bins=bin_edges, weights=np.ones(n)/n, density=False)
    counts, bins, bars=ax_hist.hist(data, bins=bin_edges, density=False, rwidth=0.85, alpha=1, label=label)
    weights=np.ones(n)/n
    errors, bin_centers=calculaerrorhisto(data, weights, bin_edges)
    
    #En realidad el output puede ser cuentas o frecuencia    
    #Esto queda feo porque pone una línea al final
    # c_KO=ax_cum.hist(intensidad_KO, bins=np.linspace(0,lim_intensidad, n_KO+1),density='true', cumulative='True', histtype="step", alpha=0.5, label='KO')
    # c_ocre=ax_cum.hist(intensidad_ocre, bins=np.linspace(0,lim_intensidad, n_ocre+1), density='true', cumulative='True', histtype="step", alpha=0.5, label='ocre')
    #Para versiones nuevas de matplotlib uso
    #c_KO=ax_cum.ecdf(intensidad_KO)
    #c_ocre=ax_cum.ecdf(intensidad_ocre)
    
    '''
    #Usaría esto si quisiese hacer una suma acumulativo con bins de 1 en 1
    freq_sin_bin_1=np.histogram(data1, bins=np.linspace(0,max_histo, n1+1), density=False, weights=np.ones(n1)/n1)
    freq_sin_bin_2=np.histogram(data2, bins=np.linspace(0,max_histo, n2+1), density=False, weights=np.ones(n2)/n2)
    cum1=np.cumsum(freq_sin_bin_1[0])
    cum2=np.cumsum(freq_sin_bin_2[0])
    '''

    #Hago el plot acumulativo
    cum=np.cumsum(freq)
    #Añado un punto al comienzo para hacer el primer step completo en el acumulativo
    cum=np.insert(cum,0,0) 
    plot_cum=ax_cum.step(bin_edges, cum, alpha=1, label=label)
    leg=ax_hist.legend()
    #ax_hist.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.02e'))
    plt.ticklabel_format(style='sci', axis='x', scilimits=(-2,3), useMathText=True)

    ylabel=ax_hist.set_ylabel('Molecule count', fontsize=14)
    ylabel=ax_cum.set_ylabel('Cumulative frequency', fontsize=14)

    return fig, ax_hist, ax_cum, bars, plot_cum, bin_centers, freq, errors, bin_edges, cum

def cambiacolorbarrashisto(bars, color, alpha=1):
    #Cambia el color de las barras de un histograma
    for bar in bars:
        bar.set_facecolor(color)
        bar.set_edgecolor(None)
        bar.set_alpha(alpha)
    #Esto lo hago para conocer las propiedades, si las necesito
    #for property, value in vars(bars[0]).items():
    #    print(property, ":", value)

def histo_x1 (data, num_bins=None, min_histo=0, max_histo=None, label='Data'):
    #Representa un histograma junto con sus errores. Sin CDF
    #Lo uso para representar, por ejemplo, los resultados de ocre, ambar o WT por separado
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    plt.rc('axes', labelsize=14)
    plt.rc('legend', fontsize=12)

    n=len(data)
    if num_bins == None:
        num_bins=np.sqrt(np.max(n))
    if max_histo != None:
        bins_histo= np.linspace(min_histo, max_histo, num_bins+1)

    fig, ax_hist = plt.subplots(1, 1)
    weights=np.ones(n)/n
    freqs, bin_edges = np.histogram(data, bins=bins_histo, weights=weights, density=False)
    counts, *_ = np.histogram(data, bins=bins_histo, density=False)
    error_counts=np.sqrt(counts)
    # bin_centers = (bin_edges[:-1] + bin_edges[1:])/2
    # bars= ax_hist.bar(bin_centers, counts, width=4, yerr=errors, label=label)
    errors_freq, bin_centers=calculaerrorhisto(data, weights, bin_edges)
    width=(bin_edges[0]-bin_edges[1])*0.85
    bars= ax_hist.bar(bin_centers, freqs, width=width, yerr=errors_freq, label=label)
    leg=ax_hist.legend()

    return fig, ax_hist, bars, bin_centers, freqs, errors_freq


def calculaerrorhisto (data, weights, bin_edges):
    #Calcula las barras de error en los histogramas
    #Los errores están en:
    #https://www.zeuthen.desy.de/~wischnew/amanda/discussion/wgterror/working.html

    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2
    errors=np.zeros(np.shape(bin_centers))

    for bin_index in range(len(bin_edges) - 1):
        # find which data points are inside this bin
        bin_left = bin_edges[bin_index]
        bin_right = bin_edges[bin_index + 1]
        in_bin = np.logical_and(bin_left <= data, data < bin_right)
        # print(f"in_bin {in_bin}")

        # filter the weights to only those inside the bin
        weights_in_bin = weights[in_bin]
        # print(f"weights_in_bin {weights_in_bin}")

        # compute the error however you want
        error = np.sqrt(np.sum(weights_in_bin ** 2))
        errors[bin_index]=error
        counts=np.sum(in_bin)
        # print(f"error {error}")

    # plot the error bars
    # plt.errorbar(bin_centers, bin_y, yerr=errors, linestyle="none")
    return errors, bin_centers
