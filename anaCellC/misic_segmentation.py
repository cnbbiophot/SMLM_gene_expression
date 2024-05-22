#MisiC. Microbe segmentation in dense colonies
#https://github.com/pswapnesh/misic/
#https://github.com/leec13/MiSiCgui
# Para analizar las imágenes de Raquel de Granada
#
# Necesito activar el entorno (en este caso lo he llamado MiSiCgui o misic)
# conda activate misic
#
# Y también instalar opencv y matplotlib
# conda install -y -c conda-forge opencv
# conda install -y -c conda-forge matplotlib
#
# Guarda una imagen segmentada de 8 bits
# Muy importante: tengo que medir la anchura de las células con imageJ (por ejemplo) e introducirlo en mean_width
# 
# jri - 30Aug19

### In case of gpu error, one might need to disabple gpu before importing MiSiC [ os.environ["CUDA_VISIBLE_DEVICES"]="-1" ]



from os import environ
from os.path import isfile, join

#Necesario para misic_segmentation
environ["CUDA_VISIBLE_DEVICES"]="-1" 
environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

from misic.misic import *
from misic.extras import *
from skimage.io import imsave,imread
from skimage.transform import resize,rescale
from skimage.filters import unsharp_mask
from skimage.exposure import adjust_gamma
import matplotlib.pyplot as plt


import cv2
import numpy as np

def segmentation(filePath, fileName, mean_width):
    # mean_width es la anchura promedio de las células en píxeles. Lo tengo que medir en imageG primero
    
    fNameOut= str() #Esta es la imagen segmentada de 8 bits
    #El  nombre tiene que ser ... loquesea_bf.tif o loquesea_fluo.tif
    #Encuentra el último _ y el punto de la extensión
    pos_us = fileName.rfind('_')
    pos_dot = fileName.rfind('.')
    fNameOut = (fileName[:pos_us]+'_segmented'+fileName[pos_us:pos_dot]+'.tif')

    # read image using your favorite package
    check_file = isfile(join(filePath, fileName))
    if check_file is False:
        print ('No se encuentra ', join(filePath, fileName))
        return fNameOut, None
    print ('Leyendo', fileName)    
    im = imread(join(filePath, fileName))
    sr,sc = im.shape

    # Parameters that need to be changed
    ## Ideally, use a single image to fine tune two parameters : mean_width and noise_variance (optional)

    #input the approximate mean width of microbe under consideration
    #Este es el valor con el que han entrenado el algoritmo
    standard_width = 9.7

    # If image is phase contrast light_background = True
    light_background = False

    # compute scaling factor
    scale = (standard_width/mean_width)

    # Initialize MiSiC
    mseg = MiSiC()

    ## preprocess using inbuit function or if you are feeling lucky use your own preprocessing
    # recomended preprocessing

    #La imagen de entrada es un uint16
    #im = adjust_gamma(im,0.25)
    #im = unsharp_mask(im,2.2,0.6)
    #im=anaCellC.laplacian_filter(im, 3)
    # la salida de adjust_gamma y unsharp mask es un float64

    if fileName[pos_us+1:pos_dot]=='bf':
        gamma=0.2
        #im = adjust_gamma(im, gamma)
        im=im.astype('float64')

    # for fluorescence images
    if fileName[pos_us+1:pos_dot]=='fluo':
        gamma=0.25
        #im = adjust_gamma(im, gamma)
        #im = unsharp_mask(im,2.2,0.6)
        im=im.astype('float64')
        im = gaussian(laplace(im),2.2)
        #im = add_noise(im,0.1)
        # OR
        # im = random_noise(im,mode = 'gaussian',var = 0.1/100.0)

    im = rescale(im,scale,preserve_range = True)

    # add local noise
    img = add_noise(im,sensitivity = 0.13, invert = light_background)

    # segment
    yp = mseg.segment(img, invert = light_background)
    yp = resize(yp,(sr,sc))
    im = resize(im,(sr,sc))

    # watershed based post processing (optional)
    #yp_post = postprocess_ws(im,yp)
    yp_post = postprocessing(im if light_background else -im, yp[:,:,0])  
    # save 8-bit segmented image and use it as you like
    imSeg=((yp_post > 0)*255).astype(np.uint8)

    print ('Saving '+fNameOut)
    imsave(join(filePath, fNameOut), imSeg)
    render_segmented_images(im, yp, yp_post, imSeg)
    return fNameOut, imSeg

def laplacian_filter(src, kernel_size):
    #imIn tiene que ser una imagen en escala de grises
    #kernel_size debe ser 2 o 3

    # [variables]
    # Declare the variables we are going to use
    #desired deptf of the destination image
    ddepth = cv2.CV_64F
          
    # [load]
    # [reduce_noise]
    # Remove noise by blurring with a Gaussian filter
    src = cv2.GaussianBlur(src, (3, 3), 0)
      
    # [laplacian]
    # Apply Laplace function
    dst = cv2.Laplacian(src, ddepth, ksize=kernel_size)
     
    # [convert]
    # converting back to uint8
    abs_dst=np.array(dst)
    #abs_dst=(2**16-1)*255*((abs_dst-np.min(abs_dst))/(np.max(abs_dst)-np.min(abs_dst)))
    #abs_dst=abs_dst.astype ('uint16')
    #Devuelve un float64
    return fNameOut, abs_dst

def render_segmented_images(imOriginal, imSeg_pre, imSeg_post, imSeg):
    #imOriginal es la imagen sin segmentar
    #imSeg_pre es la imagen segmentada tal cual
    #imSeg_post es la imagen después del watershed processing
    #imSeg es la imagen final binarizada
    '''
    fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(8, 8))
    ax = axes.ravel()
    ax[0].imshow(imOriginal, cmap=plt.cm.gray)
    ax[0].set_title('Original image')
    ax[1].imshow(imSeg_pre, cmap=plt.cm.gray)
    ax[1].set_title('Pre')
    ax[2].imshow(imSeg_post, cmap=plt.cm.gray)
    ax[2].set_title('Post')
    ax[3].imshow(imSeg, cmap=plt.cm.gray)
    ax[3].set_title('Segmented')
    '''
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(8, 4))
    ax = axes.ravel()
    ax[0].imshow(imOriginal, cmap=plt.cm.gray)
    ax[0].set_title('Original image')
    ax[1].imshow(imSeg, cmap=plt.cm.gray)
    ax[1].set_title('Segmented')
    for a in ax:
          a.axis('off')
    fig.tight_layout()
    plt.pause(0.03)
    return fig, axes
