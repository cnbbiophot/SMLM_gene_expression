#Funciones de anaCellC
# Para analizar las imágenes de Raquel de Granada
# Hay que estar en el entorno misic: 
# conda activate misic
#
#Funciones: 
# rois_to_ImageJ: Convierte los contornos en formato de lista a formato de ImageJ
# read_ImageJ_ROI
# rois_from_BF_fluo_images: Calcula los contornos a partir de las imágenes segmentadas BF y fluo
# rois_from_image: Recupera los contornos a partir de una imagen segmentada
# rois_from_image_2: 
# locaencelula: Comprueba si la localización está en una célula (o en su contorno) o no. Devuelve im_rgb e im2
# read_segmented_images: Lee las imágenes segmentadas a partir de la imagen fluo y la imagen bf
# render_segmented_image_from_ROIs: Muestra los Rois y los identifica en la imagen
# read_original_image
# renderlocalizationimages: Muestra la imagen de las células RGB y falso color
# rendercolorROIs: Pone la célula con el número de localizaciones y la representa en falso color y superpuesta a la imagen de contraste de fase
#
# jri - 15Dec23
# jri - 9May24

import matplotlib.pyplot as plt
from os.path import join

#Necesario para importar ROIs de ImageJ
#De https://pypi.org/project/roifile/
import roifile as rf

#Necesario para las demás funciones
import cv2
import numpy as np


def rois_from_BF_fluo_images (imSegB, imSegF):
    #Calcula los contornos a partir de las imágenes segmentadas BF y fluo
    # Lee la imagen segmentada y recupera los contornos
    if imSegB is None:
        imSegB = np.zeros(np.shape(imSegF), dtype='bool')
    imSegB = imSegB.astype(bool)
    if imSegF is None:
        imSegF = np.zeros(np.shape(imSegB), dtype='bool')
    imSegF = imSegF.astype(bool)
    imSeg = np.logical_or(imSegF, imSegB)
    imSeg = imSeg.astype(np.uint8) * 255
    rois=rois_from_image (imSeg)
    return rois
    
def rois_from_image (imSeg):
    # Recupera los contornos a partir de una imagen segmentada
    imSeg = imSeg.astype(np.uint8) * 255

    # Para localizar contornos. Es muy espectacular. 
    # Funciona bien para separar contornos muy cercanos, pero muchas veces los deja inconexos
    # imSeg = cv2.Canny(imSeg, 0, 200, 3)
    # print (np.shape(imSeg))
    # fig, ax = plt.subplots(1,1)
    # ax.imshow(imSeg, cmap = 'gray')
    # plt.show()

    # Contornos celulares: cnts
    #Las coordenadas de los ROIs tienen que ser int32. Si no, findContours da error
    cnts = cv2.findContours(imSeg, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_TC89_L1)
    cnts = cnts[0] if len(cnts) == 2 else cnts[1]
    numRois=len(cnts)
    #La salida de findContours es una tupla, pero yo lo convierto en lista a partir de ahora
    rois=list(cnts)
    for k in range(numRois):
        rois[k]=np.squeeze(cnts[k])

    # Find the convex hull object for each contour
    hulls = []
    for k in range(numRois):
       hull = cv2.convexHull(rois[k])
       hull= np.squeeze(hull)
       hulls.append(hull)

    # im_cnts=np.ones((np.shape(imSeg)[0], np.shape(imSeg)[1],numRois), dtype=np.uint8)*255
    # #fig, axes = plt.subplots(1,numRois)
    # #ax=axes.ravel()
    # fig, ax = plt.subplots(1,1)
    # # Dibuja los contornos y los enumera
    # for n in range (numRois):
    #     #n=4
    #     #color = (rng.randint(0,256), rng.randint(0,256), rng.randint(0,256))
    #     imtmp=np.ones((np.shape(imSeg)[0],np.shape(imSeg)[1],3), dtype=np.uint8)*255
    #     cv2.drawContours(imtmp, [hulls[n]], -1, (0, 255, 0), 1)
    #     cv2.drawContours(imtmp, [rois[n]], -1, (255, 0, 0), 1)
    #     ax.imshow(imtmp)
    #     ax.set_title(str(n))
    #     #ax.axis('off')
    #     fig.tight_layout()
    #     plt.pause(2)
    #     #plt.show()
    return hulls

def rois_from_images_2 (imSegB, imSegF):
    fig, ax = plt.subplots(1,1, figsize=(10, 10))
    ax.imshow(imSeg, cmap=plt.cm.gray)
    ax.axis('off')
    fig.tight_layout()
    plt.show()


    output=cv2.connectedComponentsWithStats(imSeg)
    labels=output[0]
    stats=output[1]
    centroids=output[2]
    print (labels)
    print (np.shape(stats))
    print (np.shape(centroids))
    return

def rois_to_ImageJ (rois, fileName):
    #Convierte los contornos en formato de lista a formato de ImageJ
    #En los contornos calculados por misic, rois es una tupla de rois
    #Convierto esos contornos en ImageJ y los guardo uno a uno
    if isinstance (rois, np.ndarray):
        roi=np.squeeze(rois)
        rIJ=rf.ImagejRoi.frompoints(roi)
        rIJ.name='01'
        rIJ.roitype=rf.ROI_TYPE.POLYGON
        rf.roiwrite(fileName, rIJ, mode='w')
    else:
        for idx, roi in enumerate(rois):
            roi=np.squeeze(roi)
            rIJ=rf.ImagejRoi.frompoints(roi)
            rIJ.name='{0:02.0f}'.format(idx)
            rIJ.roitype=rf.ROI_TYPE.POLYGON
            if idx==0:
                rf.roiwrite(fileName, rIJ, mode='w')
            else:
                rf.roiwrite(fileName, rIJ, mode='a')

def read_ImageJ_ROI (fPath, ROIFile):
    #Devuelve los rois guardados por IJ al estilo de mi programa 
    #Es decir, una lista con las coordenadas. Cada elemento de la lista es un ROI
    #Las coordenadas de los ROIs tienen que ser int32. Si no, no funciona findContours
    fName = ROIFile if ROIFile[-4:].lower() == '.zip' else ROIFile  + '.zip'
    print ('Leyendo archivo ROI: '+fName)
    rIJ=rf.roiread (join(fPath,fName))
    numRois=len(rIJ)
    rois=list()
    roiNames=list()
    for k, roi in enumerate(rIJ):
        rois.append(np.round(roi.coordinates()))
        rois[k]=rois[k].astype(np.int32)
        roiNames.append(roi.name)
    return rois, roiNames


def locaencelula(cnts, cPix, cId, cFrame, frameMax, imBF):
    # Comprueba si la localización está en una célula (o en su contorno) o no
    # cnts son los contornos
    # cPix las coordenadas de la localización
    # cId su identificación
    # cFrame contiene el número de frame en el que se produce la identificación
    # frameMax el máximo número de frames a considerar. Si frameMax es -1 analiza todos
    # imBF la imagen brightfield para la representación
    # Devuelve: numLocaCell, cumCell, cumFrame, countsPerFrame, im_rgb, im2
    # cumCell: localizaciones acumuladas por célula
    # El tamaño final de cumCell será menor que numloca porque sólo tengo en cuenta las que están en células identificadas
    # Columnas: id del contorno; id de la localización; frame de la localización; acumulativo
    # cumFrame: número de localizaciones acumulativas por frame en las células sin contar el área control
    # la primera columna es la célula (contorno); la segunda la id de la localización en Thunderstorm; 
    # la tercera el frame de la localización; la cuarta el número acumulado de localizaciones para esa célula
    # countsPerFrame: recuento de moléculas por frame para hacer un histograma estilo vida media. Lo almacena por célula
    #   La primera columna es el número de frame en el que se cuenta la molécula. Después cada columna es una célula
    #   La última columna es el área control
    # im_rgb: imagen RGB donde está marcada cada localización
    # im2: imagen en la que sumo las localizaciones para luego representarla en falso color

    numContornos=len(cnts)
    numLoca = np.shape(cPix)[0]
    numFrames = cFrame[len(cFrame) - 1]
    numLocaCell = np.zeros(numContornos, dtype=np.uint64)
    if frameMax==-1:
        frameMax=numFrames
    countsPerFrame = np.zeros((int(frameMax), numContornos+1), dtype=np.uint64) #Recuento de moléculas por frame en cada célula
    countsPerFrame[:,0]=np.arange(1, frameMax+1) #La primera columna es el número de frame en el que se cuenta la molécula
    cumCell = np.zeros((numLoca, 4), dtype=np.uint64) #localizaciones acumuladas en cada célula por frames 
    cumFrame = np.zeros((numFrames, 2), dtype=np.uint64) #localizaciones acumuladas en el total de las células por frames 
    # Número de localizaciones en una célula
    
    imSiz=np.shape(imBF)
    im_rgb = cv2.cvtColor(imBF, cv2.COLOR_GRAY2BGR)
    im2 = np.zeros(imSiz)

    l = 0  # Contador de localizaciones total en las células y el área de control
    lCell = 0  # Contador de localizaciones en las células (no área de control)
    fprev = 1  # Contador de frames
    fl = 0  # Índice del recuento de localizaciones en frames
    for idx, cp in enumerate(cPix):
        if cFrame[idx]<frameMax+1:
            # Muestra las localizaciones como puntos
            # cp contiene las coordenadas. Pueden ser maś grandes que la imagen si ésta se ha movido
            if cp[0] > imSiz[0] - 1:
                cp[0] = imSiz[0] - 2
            if cp[1] > imSiz[1] - 1:
                cp[1] = imSiz[1] - 2
            if cp[0] < 0:
                cp[0] = 0
            if cp[1] < 0:
                cp[1] = 0
            im_rgb[cp[1].astype(np.uint16), cp[0].astype(np.uint16)] = (255, 255, 0)
            im2[cp[1].astype(np.uint16), cp[0].astype(np.uint16)] = im2[cp[1].astype(np.uint16), cp[0].astype(np.uint16)] + 1
            k = -1 #Índice de contornos
            while k < numContornos - 1:
                k = k + 1
                # Busca si la coordenada (cp[1], cp[0]) está en el contorno cnts[k]
                # cnts[k] contiene el contorno k
                result = cv2.pointPolygonTest(cnts[k], cp, False)
                if result >= 0:
                    numLocaCell[k] = numLocaCell[k] + 1
                    cumCell[l, :] = [k, cId[idx], cFrame[idx], numLocaCell[k]]
                    l =l+1
                    cF=(cFrame[idx]-1).astype('uint64') #Es -1 porque los frames empiezan en 1 y los índices en countsPerFrame empiezan en 0
                    countsPerFrame[cF, k+1] = countsPerFrame[cF, k+1] + 1
                    if k < numContornos - 1:  # Si k==número de contornos es el Área de control, no la contamos como localización en célula
                        lCell = lCell+1
                    # Si ha habido un cambio de frame almacena las acumuladas en el frame anterior
                    # fprev es el frame anterior
                    # fl es el índice en cumframe donde pone las localizaciones
                    if cFrame[idx] != fprev:
                        cumFrame[fl, :] = [fprev, lCell - 1]
                        fprev = cFrame[idx]
                        fl =fl+1
                    break

    # Me quedo con las l primeras filas. Esto se llama numpy slice
    cumCell = cumCell[:l]
    # Sort by multiple keys. Super interesting
    ind = np.lexsort((cumCell[:, 3], cumCell[:, 0]))
    cumCell = cumCell[ind]
    cumFrame = cumFrame[:fl]

    im2 = 255 * ((im2 - np.min(im2)) / (np.max(im2) - np.min(im2)))
    im2 = im2.astype('uint8')
    return numLocaCell, cumCell, cumFrame, countsPerFrame, im_rgb, im2

def read_segmented_images(filePath, rootName):
    #Lee las imágenes segmentadas a partir de la imagen fluo y la imagen bf
    fileSeg_bf = rootName + '_segmented_bf.tif'
    fileSeg_fluo = rootName + '_segmented_fluo.tif'
    imSegB = cv2.imread(filePath + fileSeg_bf, cv2.IMREAD_GRAYSCALE)
    imSegF = cv2.imread(filePath + fileSeg_fluo, cv2.IMREAD_GRAYSCALE)
    return imSegB, imSegF

def render_segmented_image_from_ROIs (imOriginal, rois, roi_names):
    #imOriginal es la imagen sin segmentar
    #imSeg es la imagen final binarizada
    #El número fuera de paréntesis es el ROI y el número entre paréntesis es la célula
    imSeg=np.zeros(np.shape(imOriginal), np.uint8)
    im_bgr = cv2.cvtColor(imOriginal, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(imSeg, rois, -1, (255), -1)
    cv2.drawContours(im_bgr, rois, -1, (0, 0, 255), 1)
    #Sï quiero convertir la imagen imSeg a RGB para escribir sobre ella, etc
    #imSeg=cv2.cvtColor(imSeg,cv2.COLOR_GRAY2BGR)
    #Identifica los rois con su nombre
    for idx, roi in enumerate(rois):
        #color = (rng.randint(0,256), rng.randint(0,256), rng.randint(0,256))
        M = cv2.moments(roi)
        cx = int(M['m10'] / M['m00'])
        cy = int(M['m01'] / M['m00'])
        cv2.putText(im_bgr, roi_names[idx]+'('+str(idx)+')', (cx, cy), cv2.FONT_HERSHEY_SIMPLEX, .6, (0, 255, 0), 2)

    #Opencv organiza las imágenes en el orden BGR en vez de RGB, así que la tengo 
    #que reordenar para representar y guardar con matplotlib

    im_rgb = cv2.cvtColor(im_bgr, cv2.COLOR_BGR2RGB)
    #imSeg = cv2.cvtColor(imSeg, cv2.COLOR_BGR2RGB)
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(8, 4))
    ax = axes.ravel()
    ax[0].imshow(im_rgb)
    ax[0].set_title('Original image')
    ax[1].imshow(imSeg, cmap=plt.cm.gray)
    #ax[1].imshow(imSeg)
    ax[1].set_title('Segmented')
    for a in ax:
          a.axis('off')
    fig.tight_layout()
    plt.pause(0.03)
    return fig, axes

def read_original_image(filePath, rootName):
    # Leo la imagen brightfield. Es un uint16 y la convierto en uint8
    #Si no hay imagen BF leo la imagen fluo
    fileBF = rootName + '_bf.tif'
    fileFluo = rootName + '_fluo.tif'

    im0 = cv2.imread(join(filePath, fileBF), cv2.IMREAD_GRAYSCALE)
    if im0 is not None:
        print ('Leyendo imagen BF: ', fileBF)
        im0 = 255 * ((im0 - np.min(im0)) / (np.max(im0) - np.min(im0)))
        imOut = im0.astype('uint8')
    else:
        print ('Leyendo imagen fluo: ', fileFluo)
        im0 = cv2.imread(join(filePath, fileFluo), cv2.IMREAD_GRAYSCALE)
        imOut =np.zeros (np.shape(im0), dtype=np.uint8)
    return imOut

def renderlocalizationimages (im_bgr, im2=None):
    #Muestra la imagen de las células RGB y falso color
    # im_rgb es la imagen 
    # im2
    #im_bgr = cv2.cvtColor(imOriginal, cv2.COLOR_GRAY2BGR)
    ##cv2.drawContours(imSeg, rois, -1, (255), -1)
    #Opencv organiza las imágenes en el orden BGR en vez de RGB, así que la tengo 
    #que reordenar para representar y guardar con matplotlib
    im_rgb = cv2.cvtColor(im_bgr, cv2.COLOR_BGR2RGB)
    if im2 is None:
        fig, ax = plt.subplots(1, 1)
        ax.imshow(im_rgb)
        ax.axis('off')
        fig.tight_layout()
    else:
        fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(8, 4))
        ax = axes.ravel()
        ax[0].imshow(im_rgb)
        ax[0].set_title('RGB image')
        im=ax[1].imshow(im2, cmap=plt.cm.hot)
        ax[1].set_title('Counts')
        cbar=fig.colorbar(im, ax=ax[1])
        for a in ax:
              a.axis('off')
    fig.tight_layout()
    plt.pause(0.03)
    return fig, ax


def rendercolorROIs (im_BF, rois, nLocaCell, overlay=True, vmin=None, vmax=None, ctrlArea=True):
    #Pone la célula con el número de localizaciones para luego representarla en falso color
    # im_original es la imagen BF o fluo (sin segmentar. También sirve la segmentada, sólo es para las dimensiones)
    # imSeg es la imagen final binarizada
    # rois son los rois
    # nLocaCell es el número de localizaciones por célula
    # vmin, vmax son los valores límite de la representación de la imagen de recuento (es decir el mínimo del recuento y el máximo)
    
    im2=np.zeros(np.shape(im_BF), np.uint16)+np.min(nLocaCell)-1
    #im_bgr = cv2.cvtColor(imOriginal, cv2.COLOR_GRAY2BGR)
    #cv2.drawContours(im_original, rois, -1, (255), -1)
    
    mask = np.ones(np.shape(im2), dtype=bool) #máscara de Trues
    #Como en https://stackoverflow.com/questions/33234363/access-pixel-values-within-a-contour-boundary-using-opencv-in-python
    numRois=len(rois)
    if ctrlArea: #Si ctrlArea es true, evito usar el último roi, que es el área control
        numRois=len(rois)-1 
    for idx in range(numRois): 
        # Create a mask image that contains the contour filled in
        cv2.drawContours(im2, rois,idx, (255), thickness=-1)
        #Lo tengo que hacer en dos pasos porque quiero identificar también aquellas células que me dan 0
        #Primero relleno el contorno de 255s para poder identificarlo los puntos que lo componen
        #y luego le doy el valor del número de moléculas, aunque éste sea 0
        pts = np.where(im2 == 255) 
        im2 [pts]=(nLocaCell[idx].astype(np.float64))
        mask [pts]=False #La máscara es False si el elemento es válido

    im2_mask = np.ma.masked_array(im2, mask)
    #Sï quiero convertir la imagen imSeg a RGB para escribir sobre ella, etc
    #imSeg=cv2.cvtColor(imSeg,cv2.COLOR_GRAY2BGR)
    #Identifica los rois con su nombre

    fig, ax = plt.subplots(nrows=1, ncols=1)
    if overlay is True:
        ax.imshow(im_BF, cmap='gray', interpolation='none')
        im=ax.imshow(im2_mask, cmap='jet', alpha=0.4, interpolation='none', vmin=vmin, vmax=vmax)
        ax.axis('off')
    else:
        im=ax.imshow(im2_mask, cmap='jet', interpolation='none', vmin=vmin, vmax=vmax)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(True)

    cbar = fig.colorbar(im)
    cbar.ax.tick_params(axis='both', which='both', length=0)
    cbar.ax.set_frame_on(False)
    fig.tight_layout()
    plt.pause(0.03)
    return fig, ax
