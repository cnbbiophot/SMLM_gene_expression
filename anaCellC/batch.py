# Analiza todos los archivos que están en el documento folderlist 
# Rehace todas las figuras, etc
# Puede hacer un análisis global de todas las células 
#
# analysis: lee todos los archivos de folderlist y los analiza con el main (es decir, crea todos los archivos necesarios para luego hacer el análisis)
# readmoleculesperframe: Lee todas los archivos de recuento y devuelve el recuento por célula, el ruido, etc
# moleculesperframe: Analiza todos los archivos de foldersFile y kis ajusta
# modifyparamdat: modifica todos los param.data del folderlist
#
# O desde el intérprete o Jupyter, etc:
# import anaCellC as batch
# batch.analysis (foldersFile, 3000, frame_number) (filePath)
#
#
#
# jri - 25Jan24

import matplotlib.pyplot as plt
from os.path import join
#Necesario para las demás funciones
import numpy as np

import anaCellC.main as main
import anaCellC.utils as u
import anaCellC.plot as p
import anaCellC.IO as IO
import anaCellC.representa as r
import anaCellC.functions as f
import anaCellC.contours as contours


# importing datetime module for now()
import datetime as dt 


def analysis(foldersFile, frameFit=None, frameQuant=None):
    #Hace el main en todos los archivos de foldersFile. Eso significa leer todos los datos y generar los histogramas, etc.
    #Para poder cerrar las ventanas tengo que hacer el proceso interactivo con plt.ion()
    #Si frameFit es none, usa el valor de frameFit en param.dat
    plt.ion()
    folderList=IO.readFileList(foldersFile)
    for filePath in folderList:
        main.main (filePath, frameFit, frameQuant)
        plt.pause(1.5)
        plt.close('all')
    plt.ioff()


def readmoleculesperframe(foldersFile,frameIni=1, frameFit=None):
    #Lee todos los archivos de recuento de moléculas y devuelve las cuentas por célula, las cuentas en el frame, etc
    #Lo necesita moleculesperframe
    folderList=IO.readFileList(foldersFile)
    #Leo una vez para ver cuántos frames hay
    numFiles=len(folderList)
    f=np.zeros(numFiles, dtype=np.uint64) 
    for ndx, filePath in enumerate(folderList):
        fPath, rootName, pixSize, _, frameFitParamDat, *_ = IO.read_inputs(join(filePath,'param.dat'))
        if fPath:
            filePath=fPath
        f[ndx]=frameFitParamDat
    f0 = f[0]
    todos_iguales = (all(item == f0 for item in f))
    if todos_iguales:
        print ('Todos los frames tienen el mismo análisis:', f0)
        frameIni=int(frameIni)
        if frameFit is None:
            frameFit=int(f0)
        else:
            f0=frameFit
        print ('Analizamos hasta:', frameFit)
    else:
        print (f)
        raise Exception("ERROR: " + "No todos los frameFit son iguales")

    counts_in_frame=np.zeros((f0, 2), dtype=np.float64)
    cumLoca=np.zeros((f0, 2), dtype=np.float64)
    noise_in_frame=np.zeros(f0, dtype=np.float64)
    cumNoise=np.zeros((f0, 2), dtype=np.float64)
    #Guardo RAM para 5000 células (aunque hay menos)
    counts_in_every_cell=np.zeros((f0, 5000), dtype=np.float64)
    numCells=0
    for ndx, filePath in enumerate(folderList):
        fPath, rootName, pixSize, *_ = IO.read_inputs(join(filePath,'param.dat'))
        if fPath:
            filePath=fPath
        fileName=rootName+'_frame.xls'
        c=np.loadtxt(join(filePath,fileName), delimiter='\t', skiprows=1)
        molCell=c[frameIni-1:frameFit, 1:-1]
        counts_in_file=np.sum(molCell, axis=1) #Suma todas las localizaciones de las células
        counts_in_frame[:,1]=counts_in_frame[:,1]+counts_in_file
        numCells_in_frame=np.shape(c)[1]-2 #Es -2 Porque la primera columna de _frame.xls es el frame y la última es el área control
        #Coloca los datos de las células en counts_in_every_cell
        #counts_in_every_cell sólo tiene datos de células, no del índice de los frames
        counts_in_every_cell[:,numCells:numCells+numCells_in_frame]=molCell
        numCells=numCells+numCells_in_frame
        noise_in_file=c[frameIni-1:frameFit,-1]*numCells_in_frame
        noise_in_frame=noise_in_frame+noise_in_file #Suma las localizaciones del área control multiplicadas por el número de células en el frame

    counts_in_frame[:,0]=c[frameIni-1:frameFit,0]
    cumLoca[:,0]=c[frameIni-1:frameFit,0]
    cumLoca[:,1]=np.cumsum(counts_in_frame[:,1],axis=0)
    # Me quedo con las numCells primeras columnas. Esto se llama numpy slice
    counts_in_every_cell=counts_in_every_cell[:, :numCells]

    cumNoise[:,0]=c[frameIni-1:frameFit,0]
    cumNoise[:,1]=np.cumsum(noise_in_frame,axis=0)

    return counts_in_every_cell, counts_in_frame, noise_in_frame, cumLoca, cumNoise, numCells, f0

def moleculesperframe(foldersFile, frameIni=1, frameFit=None):
    # Analiza todos los archivos de foldersFile y kis ajusta
    print (foldersFile)
    counts_in_every_cell, counts_total, noise_total, cumLoca, cumNoise, numCells, f0=\
     readmoleculesperframe(foldersFile,frameIni, frameFit)
    if frameFit is None:
        frameFit=f0

    paramFit=u.fit_cumFrame_exp_linear (cumLoca[frameIni-1:frameFit,0], cumLoca[frameIni-1:frameFit,1])
    # paramFit=u.fit_doble_exp_complementaria (tau_KO, cumLoca[frameIni:frameFit+1,0], cumLoca[frameIni:frameFit+1,1])
    A=paramFit.x[0] #Amplitud de la exponencial acumulativa
    tau=paramFit.x[1] #Tiempo característico (en frames)
    m=paramFit.x[2] #Pendiente de la recta (número de identificaciones espúreas por frame)

    yFit=u.exp_lineal_complementaria ((A, tau, m), (cumLoca[frameIni-1:frameFit,0]).astype(np.float64))
    print('Ajuste a modelo de crecimiento acumulativo exponencial-lineal:')
    print ('Amplitud (moléculas): {0:0.2f}'.format(A))
    print ('Tiempo característico de localización (frames): {0:0.2f}'.format(tau))
    print ('Identificaciones espúreas por frame: {0:0.2f}'.format(m))
    print ('Número de células: ', numCells)
    print ('Moléculas promedio por célula: {0:0.2f}'.format (A/numCells))

    print (A, tau, m)
    print (A/numCells, tau, m/numCells)

    fig_frame= p.plotCountsPerFrame (cumLoca[frameIni-1:frameFit,:], yFit)

    plt.pause(0.03)
    plt.show()

def modifyparamdat(foldersFile, key='frameFit', value=250):
    folderList=IO.readFileList(foldersFile)
    for filePath in folderList:
        IO.modifyparamdat(filePath, key=key, value=value)
    
def valid_in_every_folder (foldersFile, valid_file, vmin=None, vmax=None, save_data=False):
    #Guarda las figuras con sólo las células válidas 
    #Para poder cerrar las ventanas tengo que hacer el proceso interactivo con plt.ion()

    plt.ion()
    folderList=IO.readFileList(foldersFile)
    valid=np.loadtxt (valid_file)
    valid=valid.astype(bool)
    molCell=np.zeros(10000, dtype=int)
    numCellsTotal=0
    for filePath in folderList:
        fPath, rootName, _, _, _, frameQuant, ROIFile, *_ = IO.read_inputs(filePath + 'param.dat')
        if fPath:
            filePath=fPath
        im0 = f.read_original_image(filePath, rootName)
        molInFile, numCellsInFile, roisInFile =r.readhist (filePath, rootName)
        rois, roiNames=f.read_ImageJ_ROI(filePath, ROIFile)
        validInFile=valid[numCellsTotal:numCellsTotal+numCellsInFile]
        print (filePath, numCellsInFile)
        print (validInFile, molInFile)
        print (numCellsInFile, np.sum(validInFile.astype(np.int16)))
        molCell[numCellsTotal:numCellsTotal+numCellsInFile]=molInFile
        numCellsTotal+=numCellsInFile
        rois=np.delete(rois, np.where(validInFile==False))
        fig_colorines, _=f.rendercolorROIs (im0, rois, molInFile[validInFile], overlay=False, vmin=vmin, vmax=vmax, ctrlArea=False)
        fig_colorines_overlay, _=f.rendercolorROIs (im0, rois, molInFile[validInFile], overlay=True, vmin=vmin, vmax=vmax, ctrlArea=False)
        fig_colorines.canvas.set_window_title('Colour rois ' + rootName)
        fig_colorines_overlay.canvas.set_window_title('Overlay '+ rootName)
        if save_data:
            #Guarda las células KO que son válidas en el análisis en un archivo para cada condición
            '''
            fName_data_out=rootName+'_valid'+'.dat'
            print('Saving', fName_data_out)
            now=dt.datetime.now()
            current_time = str(now.strftime("%d/%m/%Y, %H:%M:%S"))
            h1=fName_data_out+'\t'+current_time
            h2='Cell_number'
            h=h1+'\n'+h2
            fdata=validInFile.astype(np.int16)
            np.savetxt(join(filePath, fName_data_out), fdata, fmt='%u', header=h, delimiter='\t')
            '''
            print('Saving', rootName + '_rois_color.png')
            fig_colorines.savefig(join(filePath, rootName + '_rois_color.png'))
            print('Saving', rootName + '_rois_overlay_valid.png')
            fig_colorines_overlay.savefig(join(filePath, rootName + '_rois_overlay.png'))

        plt.pause(1.5)
        plt.close('all')
    plt.ioff()
