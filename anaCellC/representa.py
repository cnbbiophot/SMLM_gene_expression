#Funciones requeridas para las funciones de representación
# readhist: #Devuelve el total de las localizaciones en una célula en el intervalo de cuantificación
# retira_outliers
# leemolfolders #Número de moléculas en el intervalo de cuantificación en los archivos de foldersFile
# histo: representa el histograma para las figuras finales
# 
# jri Feb24
# jri 7May24

#El formato de los gráficos de:
#https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
#https://stackoverflow.com/questions/11244514/modify-tick-label-text

import anaCellC.utils as u
import anaCellC.IO as IO
import anaCellC.plot as plot

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats  
from os.path import join
# importing datetime module for now()
import datetime as dt 



def readhist (filePath, rootName):
    #Devuelve el total de las localizaciones en una célula en el intervalo de cuantificación
    fileName=rootName+'_cell_total.xls'
    dtype = np.dtype([("roi_Id", "U20"), ("cell_id", int), ("mol", int)])
    cells=np.loadtxt(filePath+fileName, dtype=dtype, delimiter='\t', skiprows=1)
    if np.shape(cells)[0]>0: 
        numCells=len(cells)-1 #Le quito el area control, que es siempre la última
        mol=np.zeros(numCells, dtype=int)
        roiNames=list()
        for k in range(numCells):
            roiNames.append(rootName+'_'+cells[k][0])
            mol[k]=cells[k][2] #El número de moléculas está en la tercera columnas 
    else:
        numCells=0
        roiNames=[]
        mol=[]
    return mol, numCells, roiNames

def retira_outliers(molCell, criterio=1.7, output=True):
    #molCell son los datos de las células KO
    #Elimino las células con número de moléculas outliers según el criterio dek iqr
    #Lo de quitar células basado en la intensidad mayor que 3 sigma solo lo debo hacer para las KO
    #porque sé que son las únicas que tienen una distribución normal

    m=np.mean(molCell)
    s=np.std(molCell)

    q75, q25 = np.percentile(molCell, [75 ,25])
    iqr = q75 - q25
    if output:
        print (m, s)
        print('Células lejos por exceso: ', np.where(molCell>m+criterio*iqr))
        print('Células lejos por defecto: ', np.where(molCell<m-criterio*iqr))
    # valores_validos=np.logical_and(molCell>(m-criterio*s), molCell<(m+criterio*s))
    valores_validos=np.logical_and(molCell>(q25-criterio*iqr), molCell<(q75+criterio*iqr))
    molCell=molCell[valores_validos]
    if output:
        print (np.mean(molCell), np.std(molCell))
    return molCell, valores_validos


def leemolfolders(foldersFile, isKO=False, criterio=1.7):
    #Número de moléculas en el intervalo de cuantificación en los archivos de foldersFile
    folderList=IO.readFileList(foldersFile)
    #Número de moléculas por célula
    molCell=np.zeros(10000, dtype=int)
    roiNames=list()
    numCellsTotal=0
    for filePath in folderList:
        fPath, rootName, _, _, _, frameQuant, *_ = IO.read_inputs(filePath + 'param.dat')
        if fPath:
            filePath=fPath
        molFile, numCellsFile, roisInFile =readhist (filePath, rootName)
        roiNames.extend(roisInFile)
        molCell[numCellsTotal:numCellsTotal+numCellsFile]=molFile
        numCellsTotal+=numCellsFile
    #Me quedo con las numCellsTotal primeras filas porque las demás son 0
    molCell=molCell[:numCellsTotal]
    if isKO == True:
        molCell, *_=retira_outliers(molCell, criterio)
    numCellsTotal=len(molCell)
    # print ('Total de células: ', numCellsTotal)
    print ('Frame de cuantificación: ', frameQuant)
    return molCell, numCellsTotal

def histo (foldersFile, sample, num_bins=None, lim_representa=None, valid_file=None, save_path=None):
    #Representa y guarda los datos para las figuras finales
    mol, num_cells= leemolfolders (foldersFile, isKO=False, criterio=1.7)
    label=sample
    #Esto hace que se grabe el texto de los svg como texto y no como path
    plt.rcParams['svg.fonttype'] = 'none'
    # frameQuant=15
    # counts_in_every_cell, num_mol, _, cum_loca, _, num_cells, *_= readmoleculesperframe(foldersFile, 1, frameQuant)
    # cum_loca_every_cell_tmp=np.cumsum(counts_in_every_cell, axis=0)
    if valid_file:
        '''
        valid=np.ones(num_cells, dtype=bool)
        x=cum_loca[:, 0]
        for f in x:
            f=f.astype(np.int16)-1
            _, v=r.retira_outliers(cum_loca_every_cell_tmp[f, :], criterio=1.7, output=False)
            valid=np.logical_and(valid, v)
        # valid=np.ones(num_cells_KO, dtype=bool) #Por si quiero ver lo que pasa si no retiro ningín outlier
        cum_loca_every_cell=cum_loca_every_cell_tmp[:, valid]
        mol=cum_loca_every_cell[-1, :]
        '''
        valid=np.loadtxt (valid_file)
        num_cells=np.sum(valid)
        valid=valid.astype(bool)
        mol=mol[valid]
    print ('Número de células: ', num_cells)
    print('Máximo número de localizaciones en células:', np.max(mol))

    if lim_representa is None:
        lim_representa=np.max(mol) #Para ocre y ambar 60
    if sample=='WT':
        lim_representa=np.max(mol) #Para WT
    # lim_representa = 5363 #Porque el WT es 5363
    print ('Límite del histograma: ', lim_representa)
    if num_bins is None:
        num_bins=int(num_cells**.5)

    #Represento 
    # fig, ax_hist0, ax_cum, *_ =p.histocumulative_x1 (mol, num_bins, 0, lim_representa, label)
    # xlabel=ax_cum.set_xlabel ('Fluorescent events per cell', fontsize=14)
    fig, ax_hist, _, bin_centers, freqs, errors_freq =plot.histo_x1 (mol, num_bins, 0, lim_representa, label)
    xlabel=ax_hist.set_xlabel ('Fluorescent events per cell', fontsize=14)

    #Análisis bootstrap
    res = stats.bootstrap((mol,), np.mean, confidence_level=0.95)
    SEM=res.standard_error
    EM_95=(res.confidence_interval.high-res.confidence_interval.low)/2

    s1="#Media, STD, SEM, 95%EM: {0:2.2f}, {1:2.2f}, {2:2.2f}, {3:2.2f}".format(np.mean(mol), np.std(mol), SEM, EM_95)
    s2="#Total de células: "+str(int(num_cells))
    s3='#Máximo número de localizaciones en células: '+ str(int(np.max(mol)))
    s4='#Lim. histograma: '+str(lim_representa)
    print (s1)

    if save_path:
        fNameOut='histo_'+label
        # print('Saving ' + join(save_path, fNameOut + 'png'))
        fig.savefig(join(save_path, fNameOut + '.png'))
        fig.savefig(join(save_path, fNameOut+'.svg'), dpi=300, format='svg')

        #Guardo los datos en un archivo de texto para Raquel:
        now=dt.datetime.now()
        current_time = str(now.strftime("%d/%m/%Y, %H:%M:%S"))

        #Guardo los datos del número de moléculas crudos
        f = open(join(save_path, 'data_' + sample + '.dat'), "w")
        f.write("#"+fNameOut+".dat"+"\t"+current_time+"\n")
        f.write(s1+'\n')        # Resultado del cálculo bootstrap en forma de comentario
        f.write(s2+'\n')
        f.write(s3+'\n')
        f.write(s4+'\n')
        f.write ('#Data\n')
        np.savetxt(f , mol, fmt='%1u', delimiter='\t')
        f.close()

        #Guardo los datos de la representación del histograma con su error
        f = open(join(save_path, fNameOut + '.dat'), "w")
        f.write("#"+'data_'+sample + ".dat\t"+current_time+"\n")
        f.write(s1+'\n')        # Resultado del cálculo bootstrap en forma de comentario
        f.write(s2+'\n')
        f.write(s3+'\n')
        f.write(s4+'\n')
        f.write ('#Bin_center\tFreq\tError\n')
        fdata=np.column_stack((bin_centers, freqs, errors_freq))
        np.savetxt(f, fdata, fmt='%0.5f', delimiter='\t')
        f.close()
    plt.show()
