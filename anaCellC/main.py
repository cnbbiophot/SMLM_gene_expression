# main de anaCellC
# Incluye la posibilidad de ejecutarse desde python con el primer argumento como nombre de archivo
# python anaCellC.py filePath
#
# O desde el intérprete o Jupyter, etc:
# import anaCellC
# anaCellC.main (filePath)
#
#O ejecutar desde anaCellC_run

# Funciones:
# main:
# segmenta: Segmenta las imágenes BF y fluo. Para ser llamada desde fuera con el filePath. 
#
# jri - 15.12.23

import anaCellC.functions as anaCellC
import anaCellC.utils as u
import anaCellC.IO as IO
import anaCellC.misic_segmentation as misic
import anaCellC.contours as contours
import matplotlib.pyplot as plt
from os.path import join

#Analiza con los siguientes parámetros:
#El nombre de la imagen de BF acaba en _bf.tif
#El nombre de la imagen de Fluorescencia acaba en _fluo.tif
#El nombre de la imagen de la secuencia acaba en _sm.tif
#Tamaño del pixel en nm. 102 para x1.6 y 4x4 binning; 81 para x1 y 2x2
#Otros parámetros:
#Ocre
#frameFit=300
#frameQuant=60
# frameMax=-1

#Ambar
#frameFit=1500
#frameQuant=200

def main (filePath, frameFit_in=None, frameQuant_in=None):
    fPath, rootName, pixSize, controlArea, frameFit, frameQuant, ROIFile, growthModel, saveData, cellWidth, colourScale = \
    IO.read_inputs(join(filePath,'param.dat'))
    if fPath:
        filePath=fPath
    if ROIFile:
        #No es necesario que el archivo de ROIs acabe en ".zip"
        print ('ROIfile: ', ROIFile)
        rois, roiNames=anaCellC.read_ImageJ_ROI(filePath, ROIFile)
        im0 = anaCellC.read_original_image(filePath, rootName)
        anaCellC.render_segmented_image_from_ROIs (im0, rois, roiNames)
    else:
        #Obtiene los ROIs de las imágenes BF y fluo
        rois_fluo, rois_BF=segmenta_2 (filePath, rootName, pixSize, cellWidth)
        rois=rois_fluo
    if colourScale is None:
        vmin=None
        vmax=None
    else:
        vmin=colourScale[0]
        vmax=colourScale[1]

    #Para cuando le doy un valor a frameFit distinto del que está en param.dat, por ejemplo frameFit=3000
    if frameFit_in is not None:
        frameFit=frameFit_in
    if frameQuant_in is not None:
        frameQuant=frameQuant_in
    contours.contours (filePath, rootName, pixSize, rois, roiNames, controlArea, frameFit, frameQuant, growthModel, vmin=vmin, vmax=vmax)
    plt.show()


def segmenta (filePath):
    # Segmenta las imágenes BF y fluo. Para ser llamada desde fuera con el filePath
    fPath, rootName, pixSize, controlArea, frameFit, frameQuant, ROIFile, growthModel, plotFits, cellWidth= \
    IO.read_inputs(join(filePath, 'param.dat'))
    if fPath:
        filePath=fPath
    _,_=segmenta_2 (filePath, rootName, pixSize, cellWidth)

def segmenta_2 (filePath, rootName, pixSize, cellWidth):
    #Segmentación. Puede hacerlo sobre las imágenes _fluo o las imágenes _bf
    #Hago la segmentación de las imágenes
    #pixSize es en nm
    #cellWidth es la anchura celular en um
    # mean_width es la anchura promedio de las células en píxeles es su anchura medida por Raquel en nm entre el tamaño del píxel
    # Lo tengo que medir en ImageG primero
    mean_width = 1000*cellWidth/pixSize 
    fileName=rootName+'_bf.tif'
    fNameOut, imSegB = misic.segmentation(filePath, fileName, mean_width)
    fileName=rootName+'_fluo.tif'
    fNameOut, imSegF = misic.segmentation(filePath, fileName, mean_width)
    rois_fluo=[]
    rois_BF=[]
    #imSegB, imSegF = anaCellC.read_segmented_images(filePath, rootName)
    if imSegF is not None:
        rois_fluo = anaCellC.rois_from_image (imSegF) 
        ROIFile=rootName+'_ROI_fluo.zip'
        print ('Guardando los rois en el archivo '+ROIFile)
        anaCellC.rois_to_ImageJ (rois_fluo, join(filePath,ROIFile))
    if imSegB is not None:
        rois_BF = anaCellC.rois_from_image (imSegB)
        ROIFile=rootName+'_ROI_BF.zip'
        print ('Guardando los rois en el archivo '+ROIFile)
        anaCellC.rois_to_ImageJ (rois_BF, join(filePath,ROIFile))
    print ('OK')
    return rois_fluo, rois_BF

if __name__ == '__main__':
    import sys
    main (sys.argv[1])

#main (filePath)