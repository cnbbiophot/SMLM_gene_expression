# Lee los datos del eje mayor y menor de las ROIS "ajustadas" de Raquel
# Así calculo el volumen de una bacteria como una geometría esferocilíndrica
#
#
# Basado en copia_archivos_ROI y anaCellC_utils
# jri - 28.11.23
# jri - 19.5.24

from os import listdir
from os.path import isdir, isfile, join
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats  



filePath='/mnt/data/jri/lab/!Experimental/2023/mEos2/ROIs_bacteria_size/'

file_list = [f for f in listdir(filePath) if isfile(join(filePath, f))]
dir_list = [f for f in listdir(filePath) if isdir(join(filePath, f))]
delimiter=','

def calculaerror (m):
    res = stats.bootstrap((m,), np.mean, confidence_level=0.95)
    SEM=res.standard_error
    EM_95=(res.confidence_interval.high-res.confidence_interval.low)/2
    return SEM, EM_95

majorAx=np.zeros(10000, dtype='float64')
minorAx=np.zeros(10000, dtype='float64')
k=0
for f in file_list:
    fileName=join(filePath, f)
    #print (fileName)
    fid = open(fileName, "r")
    l=0
    for line in fid:
        if l>0:
            if line[0] == "#" or len(line) < 2:
                continue
            cleanLine = line.strip().split(delimiter)
            majorAx[k]=float(cleanLine[2])
            minorAx[k]=float(cleanLine[3])
            k=k+1
        l=l+1
    fid.close()

# Me quedo con las k primeras filas. Esto se llama numpy slice
majorAx = majorAx[:k]
minorAx = minorAx[:k]
V=np.pi*minorAx**2*(majorAx-minorAx/3)/4

#Hago el bootstrapping
SEM_majorAx, EM_95_majorAx=calculaerror(majorAx)
SEM_minorAx, EM_95_minorAx=calculaerror(minorAx)
SEM_V, EM_95_V=calculaerror(V)

print ('Número de células: ', k)
print ('L (um), SEM, 95%SEM: ', np.mean(majorAx), np.std(majorAx)/k**.5, SEM_majorAx, EM_95_majorAx)
print ('w (um), SEM, 95%SEM: ', np.mean(minorAx), np.std(minorAx)/k**.5, SEM_minorAx, EM_95_minorAx)
print ('Ratio, SEM', np.mean(majorAx/minorAx), np.std(majorAx/minorAx)/k**.5)
print ('Volumen (fL), SEM, 95%SEM: ', np.mean(V), np.std(V)/k**.5, SEM_V, EM_95_V) 

fig_hist, ax_hist = plt.subplots(1, 1)
output=ax_hist.hist(V, bins=20, density='true', align='mid', rwidth=0.85)
plt.show()
