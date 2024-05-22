# Traza los contornos de la imagen de las bacterias binarizada usando opencv
# Para saber si las coordenadas de una molécula están dentro de un contorno uso:
# https://stackoverflow.com/questions/13786088/determine-if-a-point-is-inside-or-outside-of-a-shape-with-opencv
# https://stackoverflow.com/questions/68047866/efficient-way-to-find-coordinates-of-connected-blobs-in-binary-image
#
# Coge las imágenes segmentadas de BF y las segmentadas de fluorescencia y les hace un or
# Por ahora hay que indicar a mano las coordenadas de la célula control, que es la última de todas
#
# Más documentación:
# https://docs.opencv.org/3.4/dd/d49/tutorial_py_contour_features.html
#
# Funciones:
# contours
# reorganiza_numLocaCell: Añade una columna a nLocaCell con el índice del contorno
# dibujaContornos: Dibuja los contornos indicando el número de moléculas identificadas
# cum_exp_linear_growth
# cum_linear_growth
# fitGrowthModel: Ajuste del acumulativo en el frame

#Incluye ejemplos de cómo grabar texto y datos a la vez

#
# jri - 30Aug23
# jri - 10Jan24 - Añado el ajuste acumulativo por partes (recta primero, exponencial acumulativa después)
# jri - 25Jan24 - Añado el histograma de recuento de moléculas por frame al estilo exponential decay
# jri - 15Feb24 - Guarda imágenes con colores superpuestos




# Standard imports
import cv2
import numpy as np
import matplotlib.pyplot as plt
from os.path import join
#import random as rng

import anaCellC.functions as f
import anaCellC.utils as u
import anaCellC.plot as p
import anaCellC.IO as IO

def contours (
    filePath, rootName, pixSize, rois, roiNames, ctrlArea, frameFit, frameQuant, cumGrowthModel, save_data=True, vmin=None, vmax=None, numDetections=1):
    #frameFit es el máximo frame para hacer el ajuste
    #frameQuant es el máximo frame para hacer la cuantificación teniendo en cuenta el tiempo característico
    # fileLoca = rootName+'.xls' o +'.csv'
    #numDetections es el número de frames en los que se detecta una molécula
    fNameOut = rootName + '_rgb.tif'
    fNameCellCum = rootName + '_cell_cum.xls' #Acumulativa por frame y célula
    fNameCellTotal = rootName + '_cell_total.xls' #Total por célula en el tiempo de cuantificación
    fNameFrameCum = rootName + '_frame_cum.xls' #Acumulativa de todas las cálulas por frame. Incluye el ajuste
    fNameCountPerFrame = rootName + '_frame.xls' #Total moléculas por célula y frame

    im0 = f.read_original_image(filePath, rootName)
    # Lee las coordenadas de las localizaciones en pixeles (cPix), 
    #su identificación (cId) y el frame en el que se encuentran
    #cId, cFrame, cPix, *_ = u.read_loca_coordinates(filePath, rootName, pixSize)
    cId, cFrame, cPix, _, _, _, _, _, _, cDetections = IO.read_loca_coordinates(filePath, rootName, pixSize)
    loca_mas_de_una=cDetections>=numDetections
    cId=cId[loca_mas_de_una]
    cFrame=cFrame[loca_mas_de_una]
    cPix=cPix[loca_mas_de_una]
    numLoca = np.shape(cPix)[0]
    numFrames = cFrame[len(cFrame) - 1]
    numCells=len(rois) #Número de células
    if frameFit==-1:
        frameFit=numFrames
    print ('Número de localizaciones: ', numLoca)
    print ('Último frame: ', cFrame[numLoca-1])
    print('Número de frames: ', numFrames)

    # Contornos de área control: rois_ctrl
    # Creo un rectángulo que me sirva de área control y la añado a las células (contornos) detectados
    x0 = ctrlArea[0]; y0 = ctrlArea[1]; w = ctrlArea[2]; h = ctrlArea[3]
    cnt_ctrl = np.array([[[x0, y0]], [[x0 + w, y0]], [[x0 + w, y0 + h]], [[x0, y0 + h]]], dtype=np.int32)
    # Las áreas control van siempre a partir de la última célula detectada
    rois.append(cnt_ctrl)
    numContornos = len(rois) #El último contorno es el área de control

    #Identifica la célula a la que corresponde cada localización 
    #Para el número de frames de ajuste (_f) y de  cuantificación (_q)
    nLocaCell_f, cumCell_f, cumFrame_f, countsPerFrame, im_bgr_f, im2_f= f.locaencelula (rois, cPix, cId, cFrame, frameFit, im0)
    nLocaCell_q, cumCell_q, cumFrame_q, _, im_bgr_q, im2_q= f.locaencelula (rois, cPix, cId, cFrame, frameQuant, im0)
    # cumFrame da el número de localizaciones acumulativas por frame en las células sin contar el área control
    # la primera columna es la célula (contorno); la segunda la id de la localización en Thunderstorm; 
    # la tercera el frame de la localización; la cuarta el número acumulado de localizaciones para esa célula

    #nLocaCell es el total localizaciones en cada célula. La última célula es el área control
    #En nLocaCell el índice es el número de contorno. Lo hago explícito numLocaCell:
    #La primera columna es el número de contorno. La segunda el número de localizaciones por célula
    numLocaCell_f=reorganiza_numLocaCell (numContornos, nLocaCell_f)
    numLocaCell_q=reorganiza_numLocaCell (numContornos, nLocaCell_q)

    # Dibuja los contornos y los enumera
    #Por ahora sólo con la imagen de cuantificación 
    im_bgr, im2 = dibujaContornos (rois, numLocaCell_q, im_bgr_q, im2_q)
    fig_rgb, _=f.renderlocalizationimages (im_bgr) #La figura de las localizaciones azules
    # fig_rgb, _=f.renderlocalizationimages (im_bgr, im2)
    fig_colorines, _=f.rendercolorROIs (im0, rois, nLocaCell_q, overlay=False, vmin=vmin, vmax=vmax) #Los rois de colores
    fig_colorines_overlay, _=f.rendercolorROIs (im0, rois, nLocaCell_q, overlay=True, vmin=vmin, vmax=vmax) #Los rois de colores superpuestos a la imagen de contraste de fase
    fig_rgb.canvas.set_window_title(rootName)
    fig_colorines.canvas.set_window_title('Colour rois ' + rootName)
    fig_colorines_overlay.canvas.set_window_title('Overlay '+ rootName)
    plt.pause(0.03)

    #Hago el ajuste del acumulativo en el frame:
    paramFit, yFit, fit = fitGrowthModel(cumGrowthModel, cumFrame_f, numCells, frameFit, numLocaCell_f, frameQuant, numLocaCell_q)

    fig_cum = p.plot2Graphs (numContornos, cumCell_f, cumFrame_f, yFit, frameFit, frameQuant) # Acumulativo de células y frame
    fig_frame = p.plotcumframe  (numContornos, cumFrame_f, yFit, frameFit, frameQuant) #Acumulativo en el frame (sólo esa figura)
    fig_cum.canvas.set_window_title('Cumulative cell - ' + rootName)
    fig_frame.canvas.set_window_title('Cumulative frame - ' + rootName)
    # fig_hist = p.plotHisto (numLocaCell_q) #Histograma de frame
    # fig_TR_events, *_ = p.plotCountsPerFrame (countsPerFrame) #Detección resuelta en el tiempo
    # fig_hist.canvas.set_window_title(rootName)
    # fig_TR_events.canvas.set_window_title('Time resolved events - ' +  rootName)

    if save_data:
        #Guardo todo:
        #Recuento acumulativo de localizaciones en células. La última es el área control
        print('Saving', fNameCellCum)
        h = 'cell_id\tloca_id\tframe\tacum'
        np.savetxt(join(filePath, fNameCellCum), cumCell_q, header=h,fmt='%1u', delimiter='\t')
        #Recuento acumulativo de la suma de localizaciones en células en el frame. Incluye el ajuste
        print('Saving', fNameFrameCum)
        h = 'frame\tacum\tfit'
        cumFrame_f_fit=np.column_stack((cumFrame_f,  yFit))
        np.savetxt(join(filePath, fNameFrameCum), cumFrame_f_fit, header=h, fmt='%1u', delimiter='\t')
        #Total de localizaciones en células en el frame en el intervalo de cuantificación (q) con roi_id
        print('Saving', fNameCellTotal)
        h = 'roi_id\tcell_id\ttotal'
        roiNames.append('control')
        roi_q =  np.column_stack((roiNames,  numLocaCell_q))
        np.savetxt(join(filePath, fNameCellTotal), roi_q, header=h, fmt='%s', delimiter='\t')
        #Acumulativo de localizaciones en cada célula incluida el área control
        print('Saving', fNameCountPerFrame)
        h = 'frame'
        for k in range(numContornos-1):
            h=h+'\tcell_'+('{0:02.0f}').format(k)
        h=h+'\tctrl_area'
        np.savetxt(join(filePath, fNameCountPerFrame), countsPerFrame , header=h, fmt='%1u', delimiter='\t')

        #Guardo las figuras
        #Imagen BGR
        print('Saving', fNameOut)
        cv2.imwrite(join(filePath, fNameOut), im_bgr)
        print('Saving', rootName + '_cum.png')
        fig_cum.savefig(join(filePath, rootName + '_cum.png'))
        print('Saving', rootName + '_frame.png')
        fig_frame.savefig(join(filePath, rootName + '_frame.png'))
        print('Saving', rootName + '_rois_color.png')
        fig_colorines.savefig(join(filePath, rootName + '_rois_color.png'))
        print('Saving', rootName + '_rois_color_2.png')
        fig_colorines_overlay.savefig(join(filePath, rootName + '_rois_overlay.png'))
        # print('Saving', rootName + '_hist.png')
        # fig_hist.savefig(join(filePath, rootName + '_hist.png'))
        # print('Saving', rootName + '_TRframe.png')
        # fig_TR_events.savefig(join(filePath, rootName + '_TRevents.png'))

    return paramFit


def reorganiza_numLocaCell (numContornos, nLocaCell):
    #Añade una columna a nLocaCell con el índice del contorno
    ntmp = np.zeros((numContornos, 2), dtype=np.uint64)
    ntmp[:, 0] = np.arange(numContornos)
    ntmp[:, 1] = nLocaCell
    numLocaCell = ntmp
    return numLocaCell

def dibujaContornos (rois, numLocaCell, im_bgr, im2):
    #Dibuja los contornos indicando el número de moléculas identificadas
    for idx, c in enumerate(rois):
        #color = (rng.randint(0,256), rng.randint(0,256), rng.randint(0,256))
        cv2.drawContours(im_bgr, [c], -1, (0, 0, 255), 1)
        cv2.drawContours(im2, [c], -1, (0, 0, 255), 1)
        #print (idx, numLocaCell[idx,1])
        M = cv2.moments(c)
        cx = int(M['m10'] / M['m00'])
        cy = int(M['m01'] / M['m00'])
        #cv2.putText(im_bgr, str(numLocaCell[idx, 1])+'('+str(idx)+')', (cx+12, cy), cv2.FONT_HERSHEY_SIMPLEX, .6, (0, 255, 0), 2)
        cv2.putText(im_bgr, str(numLocaCell[idx, 1]), (cx + 12, cy), cv2.FONT_HERSHEY_SIMPLEX, .6, (0, 255, 0), 2)
        cv2.putText(im2, str(numLocaCell[idx, 1]), (cx + 12, cy), cv2.FONT_HERSHEY_SIMPLEX, .6, (0, 255, 0), 2)
    #cv2.circle(im_rgb, point, 8, (100, 100, 255), -1)
    return im_bgr, im2

def cum_exp_linear_growth(cumFrame, numCells, frameFit, numLocaCell_f, frameQuant, numLocaCell_q):
    #frameNum es el índice del frame; cumLoca es el número de localizaciones en ese frame
    #frameFit es el número de frames que ajusta; frameQuant el número de frames de cuantificación
    #numLocaCell la matriz de localizaciones p en el intervalo de ajuste (_f) y de cuantificación (_q) 
    frameNum=cumFrame[:, 0]
    cumLoca=cumFrame[:, 1]
    totalLocaCell_f=np.sum(numLocaCell_f[:,1])
    totalLocaCell_q=np.sum(numLocaCell_q[:,1])
    #Lo siguiente es si ajusto con todos los parámetros libres
    paramFit=u.fit_cumFrame_exp_linear (frameNum, cumLoca)
    A=paramFit.x[0] #Amplitud de la exponencial acumulativa
    tau=paramFit.x[1] #Tiempo característico (en frames)
    m=paramFit.x[2] #Pendiente de la recta (número de identificaciones espúreas por frame)
    
    yFit=u.exp_lineal_complementaria ((A, tau, m), frameNum.astype(np.float64))
    print('Ajuste a modelo de crecimiento acumulativo exponencial-lineal:')
    print ('Amplitud (moléculas): {0:0.2f}'.format(A))
    print ('Tiempo característico de localización (frames): {0:0.2f}'.format(tau))
    print ('Identificaciones espúreas por frame: {0:0.2f}'.format(m))

    #Resumen    
    print ('Número de células: ', numCells)
    print ('')
    print ('Intervalo de ajuste --------------------------------------')
    print ('Número de frames analizado: ', frameFit)
    print('Total localizaciones en células: ', totalLocaCell_f)
    print('Número localizaciones promedio por célula: {0:0.2f}'.format(totalLocaCell_f/numCells))
    print ('Total indentificaciones espúreas: {0:0.0f}'.format (m*frameFit))
    print ('Moléculas promedio por célula: {0:0.2f}'.format (A/numCells))
    print ('Total indentificaciones espúreas por célula: {0:0.2f}'.format (m*frameFit/numCells))

    print ('Intervalo de cuantificación ------------------------------')
    print ('Número de frames analizado: ', frameQuant)
    print('Total localizaciones en células: ', totalLocaCell_q)
    print('Número localizaciones promedio por célula: {0:0.2f}'.format(totalLocaCell_q/numCells))
    print ('Total indentificaciones espúreas: {0:0.0f}'.format (m*frameQuant))
    print ('Total indentificaciones espúreas por célula: {0:0.2f}'.format (m*frameQuant/numCells))
    return paramFit, yFit

def cum_linear_growth(cumFrame, numCells, frameFit, numLocaCell):
    #frame es el índice del frame; cumLoca es el número de localizaciones en ese frame
    #frameFit es el número de frames que ajusta; numLocaCell la matriz de localizaciones
    frameNum=cumFrame[:, 0]
    cumLoca=cumFrame[:, 1]
    totalLocaCell=np.sum(numLocaCell[:,1])
    paramFit=u.fit_cumFrame_linear (frameNum, cumLoca)
    yFit=u.straightline (paramFit.x, frameNum.astype(np.float64))
    print('Ajuste a modelo de crecimiento acumulativo lineal:')
    print ('Identificaciones por frame (m): {0:0.2f}'.format(paramFit.x[0]))
    print ('Identificaciones frame inicial (y0): {0:0.2f}'.format(paramFit.x[1]))
    print ('Identificaciones experimentales frame inicial: {0:0.0f}'.format(cumLoca[0]))

    #Resumen    
    print ('Número de células: ', numCells)
    print ('')
    print ('Intervalo de ajuste --------------------------------------')
    print ('Número de frames analizado: ', frameFit)
    print('Total localizaciones en células: ', totalLocaCell)
    print('Número localizaciones promedio por célula: {0:0.2f}'.format(totalLocaCell/numCells))
    return paramFit, yFit

def fitGrowthModel(cumGrowthModel, cumFrame_f, numCells, frameFit, numLocaCell_f, frameQuant, numLocaCell_q):
    #Hago el ajuste del acumulativo en el frame:
    if cumGrowthModel.casefold()=='exp-linear':
        fit, yFit=cum_exp_linear_growth (cumFrame_f, numCells, frameFit, numLocaCell_f, frameQuant, numLocaCell_q)
        paramFit=fit.x
    elif cumGrowthModel.casefold()=='linear':
        fit, yFit=cum_linear_growth (cumFrame_f, numCells, frameFit, numLocaCell_f)
        paramFit=fit.x
    else:
        print ("No hay modelo de ajuste")
        fit=[]
        paramFit=[]
        yFit=[]
    return paramFit, yFit, fit
