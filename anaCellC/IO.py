#Funciones de lectura que usa anaCellC
#Pueden ser usadas en otros programas
#Funciones:
# readInputFile: lee las entradas de un archivo (genérico)
# read_inputs: lee los parámetros de un arvhivo (específicos de anaCellC)
# read_loca_coordinates: lee las coordenadas de un archivo Thunderstorm 
# readFileList: lee un archivo que contiene los nombres de los archivos para analizar conjuntamente (lo uso para los histogramas, etc)
# nonBlankLines: usada por readFileList
# modify: modifica un archivo de texto (para modificar param.dat)
# modifyparamdat: modifica un parámetro en param.dat
#
# jri - 31.1.24
# jri - He añadido los inputs como global en el archivo


from os.path import isfile, join
import numpy as np
from scipy import optimize
import re #Para modificar el archivo param.dat

# importing datetime module for now()
import datetime as dt 

required_inputs=["rootName", "pixelSize", "controlArea", "frameFit", "frameQuant"]
optional_inputs=["ROIFile", "saveData", "filePath", "growthModel", "cellWidth", "colourScale"]

def readInputFile(filename, inputs=[], optionalInputs=[], default=[], delimiter=" ", ignore="#"):
    file = open(filename, "r")
    allInputsDict = {}
    for line in file:
        if line[0] == "#" or len(line) < 2:
            continue
        cleanLine = line.strip().split(delimiter)
        #allInputsDict[cleanLine[0]] = float(cleanLine[1])
        allInputsDict[cleanLine[0]] = cleanLine[1]
    allInputsNames = list(allInputsDict.keys())
    # Delete not expected inputs
    for name in allInputsNames:
        if name not in inputs and name not in optionalInputs:
            del allInputsDict[name]
    # Check if all the expected inputs are in the file
    for inputName in inputs:
        if inputName not in allInputsDict.keys():
            raise Exception("ERROR: " + inputName + " not found in " + filename)
    if (len(optionalInputs) == len(default)):
        for optInp, defaultValue in zip(optionalInputs, default):
            if optInp not in allInputsNames:
                allInputsDict[optInp] = defaultValue
    return allInputsDict

def read_inputs(fileName  = 'param.dat'):
    inputs=readInputFile(fileName, inputs = required_inputs, optionalInputs = optional_inputs, delimiter = ":")
    #inputs es un diccionario
    rootName = inputs["rootName"]
    rootName=rootName.strip(" '")
    pixelSize = float(inputs["pixelSize"])
    cA= inputs["controlArea"]
    cA=cA.strip(" []")
    controlArea=np.fromstring(cA, dtype=int, sep=',')
    frameFit= float(inputs["frameFit"])
    frameQuant= float(inputs["frameQuant"])
    ROIFile = str() #Cadena vacía
    saveData = True
    filePath = str() #Cadena vacía
    if "ROIFile" in inputs.keys():
        ROIFile = inputs["ROIFile"]
        ROIFile = ROIFile.strip(" '")
    if "saveData" in inputs.keys():
        saveData = bool(inputs["saveData"])
    if "filePath" in inputs.keys():
        filePath = inputs["filePath"]
        filePath = filePath.strip(" '")
    growthModel='exp_linear'
    if "growthModel" in inputs.keys():
        growthModel = inputs["growthModel"]
        growthModel = growthModel.strip(" '")
    cellWidth=0.967
    if "cellWidth" in inputs.keys():
        cellWidth = float(inputs["cellWidth"])
    colourScale=None
    if "colourScale" in inputs.keys():
        cS= inputs["colourScale"]
        cS=cS.strip(" []")
        colourScale = np.fromstring(cS, dtype=int, sep=',')


    return filePath, rootName, pixelSize, controlArea, frameFit, frameQuant, ROIFile, growthModel, saveData, cellWidth, colourScale

def readFileList (foldersFile):
    #Lee la lista de archivos en un archivo con paths
    fid = open(foldersFile, "r")
    # reading the file 
    data = fid.read()
    # trimmin spaces and ' and replacing end splitting the text when newline ('\n') is seen. 
    f = data.split("\n") 
    folderList = nonBlankLines(f)
    fid.close() 
    return folderList

def nonBlankLines(list_in):
    #usada por readfilelist
    list_out=[]
    for l in list_in:
        line = l.rstrip(" '")
        if line:
            if line[0]!="#":
                list_out.append(l)
    return list_out


def read_loca_coordinates(filePath, rootName, pixSize):
    #Lee las coordenadas de un archivo ThunderStorm
    #pixSize en nm. Devuelve el frame de la localización en la columna 2 
    #[Constants]
    COLID=0
    COLFRAME=1
    COLX=2
    COLY=3
    COLSIGMA=4
    COLINTENSITY=5
    COLOFFSET=6
    COLBKGSTD=7
    COLCHI2=8
    COLUNCERTAINTY=9
    COLDETECTIONS=10

    fileName=rootName+'.xls'
    check_file = isfile(join(filePath,fileName))
    if check_file is True:
        locData=np.loadtxt(join(filePath,fileName), delimiter='\t', skiprows=1)
    else:
        fileName=rootName+'.csv'
        locData=np.loadtxt(join(filePath,fileName), delimiter=',', skiprows=1)

    numLoc=np.shape(locData)[0]

    #coordenadas de las localizaciones en pixeles y el frame
    cPix=np.empty([numLoc, 2])
    cPix[:, 0]=locData[:, COLX]
    cPix[:, 1]=locData[:, COLY]
    cPix=np.floor(cPix/pixSize)
    cId=locData[:, COLID]
    cId=cId.astype(np.uint64)
    cFrame=locData[:, COLFRAME]
    cFrame=cFrame.astype(np.uint64)
    cSigma=locData[:, COLSIGMA] #Sigma de la PSF en nm
    cIntensity=locData[:, COLINTENSITY] #Intensidad integrada en fotones
    cOffset=locData[:, COLOFFSET] #Intensidad offset en fotones
    cBkgStd=locData[:, COLBKGSTD] 
    cChi2=locData[:, COLCHI2]
    cUncertainty=locData[:, COLUNCERTAINTY] #Incertidumbre de la localización
    cDetections=locData[:, COLDETECTIONS] #Número de detecciones de la misma molécula
    cDetections=cDetections.astype(np.uint64)
    # print (np.type(cId))
    # print (np.type(cPix))
    # print (np.type(cDetections))
    return cId, cFrame, cPix, cSigma, cIntensity, cOffset, cBkgStd, cChi2, cUncertainty, cDetections



def modify(filepath, from_, to_):
#Modifica param.dat
# De https://www.askpython.com/python/built-in-methods/modify-text-file-python
    file = open(filepath,"r+")
    text = file.read()
    pattern = from_
    splitted_text = re.split(pattern,text)
    modified_text = to_.join(splitted_text)
    with open(filepath, 'w') as file:
        file.write(modified_text)

def modifyparamdat(filePath, key, value):
    #Modifica param.dat
    fileName=join(filePath, 'param.dat')
    print (fileName)
    inputs=readInputFile(fileName, inputs = required_inputs,
                             optionalInputs = optional_inputs, delimiter = ":")
    #Comprueba si value es un número o un string
    if isinstance(value, str) == False:
        value=str(value)
    if key in inputs:
        print (key, inputs[key])
    if not value:     #Si value está vacía y el archivo contiene key, la borra
        if key in inputs:
            del inputs[key]
            print (key + ' deleted')
    else: #Si value no está vacía le cambia el valor
        inputs[key]=' '+value
        print (key, inputs[key])
    h=list()
    h.append('#Parámetros obligatorios')
    h.append('#El separador es un espacio')
    h.append('#Los parámeteros acaban en ":"')
    h.append('#rootName: nombre raíz de todos los nombres de los archivos')
    h.append('#pixelSize en nm (102 para binning 4 y x1.6, 81 para binning 2)') 
    h.append('#controlArea (en pixeles): Formato: x0,y0,width,height')
    h.append('#Valores típicos: w=7,h=25 para binning 4x1.6; w=11,h=40 para binning 2;') 
    h.append('#frameFit: último frame para hacer el fitting. Si es -1 ajusta todos')
    h.append('#frameQuant: último frame para hacer la cuantificación')
    h.append('#Parámetros opcionales: savePlots y filePath, ROIFile, growthModel, cellWidth')
    h.append('#growthModel es el modelo de crecimiento acumulativo puede ser exp_linear, linear (u otros). Por defecto es exp_linear')
    h.append('#cellWidth es la longitud de la célula en um. Si no se especifica es el que midió Raquel: 0.967 um')
    h.append('#colourScale son los límites de la escala de colores en las imágenes de recuento (falso color y superpuesta).')
    h.append('#Si no se especifica colourScale pone los límites mínimo y máximo.')

    now=dt.datetime.now()
    current_time = str(now.strftime("%d/%m/%Y, %H:%M:%S"))

    f = open(fileName,"w+")
    f.write("#"+current_time+"\n")
    for s in h:
        f.write(s+'\n')
    f.write('\n')
    for k,v in inputs.items():
        f.write (k+':'+v+'\n')
    f.close()



