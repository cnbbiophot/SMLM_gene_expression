# Instrucciones para anaCellC
#
# jri - 20.12.23
# jri - 7.5.24

# Instrucciones para instalar misic
# Crea el entorno para misic
conda create -y --name misic git python=3.7
# Instala misic
conda activate misic
pip install misic
conda install -y -c conda-forge opencv
conda install -y -c conda-forge matplotlib
pip install -U roifile[all]

# Usar Jupyter en el entorno nuevo:
conda activate base
(base) conda install -c conda-forge nb_conda_kernels
(base) conda activate misic
(misic) conda install ipykernel
(misic) conda deactivate
# Y luego ejecuto
(base) jupyter-notebook

#Si no funciona:
https://github.com/jupyter/notebook/issues/5014
pip install --upgrade jupyter_client


# Si alguna de las librerías se ha quedado desfasada
pip show matplotlib
pip install --upgrade matplotlib 



#######
# Instrucciones para utilizar anaCellC
1. Abres Jupyter-notebook desde el entorno base (conda activate base) y escoges el entorno misic
2. Cambias a la carpeta en la que hayas guardado anaCellConcentration. Por ejemplo:
cd /my_python/anaCellConcentration/ (shift+Enter)
3. import anaCellC (shift+Enter)
3b. Si da error de hfd5, escribir en en entorno misic: conda install -y -c conda-forge hdf5=1.10.5
4. Escribes: filePath='carpeta en la que están los archivos para analizar' Por ejemplo:
filePath='/mnt/data/jri/lab/!Experimental/2023/mEos2/Oct23/Ocre_epi/14_231124_01c/'
5. En esa carpeta tiene que haber un archivo param.dat. Tienes el modelo en anaCellConcentration. 
6. Rellenas el archivo. rootname es el nombre raíz del experimento. En ROIFile pones ROIFile: rootname_ROi.zip, donde rootname es el nombre raíz del experimento
7. Haces la segmentación automática: anaCellC.segmenta (filePath) (shift+Enter)
8. Se crean dos archivos rootname_roi_bf.zip y rootname_roi_fluo.tif
9. Corriges ROI Manager de ImageJ y lo guardas como rootname_ROI.zip
11. Haces la cuantificación: anaCellC.main (filePath) (shift+Enter)
12. Guardas el labbook y sigues trabajando. 

# Instrucciones para trabajar en lotes
Operaciones previas:
1. Descargar el folderlist de cada condición
2. Corregir para Windows si es necesario
3. Descargar el valid de las KO

Desde Jupyter-notebook, entorno misic
1. Cambia a la carpeta en la que esté anaCellConcentration: 
cd /my_python/anaCellConcentration/ (shift+Enter)
2. import anaCellC.batch
3. import anaCellC.representa
4. Crea una carpeta para cada condición de cuantificación
5. Copia las carpetas con los archivos que vayas a analizar
6. Crea un foldersFile en el que indicas los archivos que se analizarán
7. batch.modifyparamdat(foldersFile, key='frameFit', value=6000) #Para WT
7. batch.modifyparamdat(foldersFile, key='frameFit', value=2800) #Para Ocre y KO
8. batch.modifyparamdat(foldersFile, key='frameQuant', value=el que quieras)
9  batch.modifyparamdat(filePath, key='colourScale', value='min, max')
9.  Para borrar la palabra clave usas: batch.modifyparamdat(filePath, key='colourScale', value='')
10. anaCellC.batch.analysis(foldersFile)
11. Para representar sólo las vcélulas válidas en el caso de las KO:
12. batch.valid_in_every_folder (foldersFile, valid_file, vmin=0, vmax=20, save_data=True)
13. Para representar el histograma: save_path es el path a la carpeta donde está folderlist. Si no se pone no graba
13. representa.histo (foldersFile, sample, valid_file=valid_file, num_bins=12, lim_representa=lim_representa, save_path=save_path)
14. 
14. Más instrucciones anacell_batch_run.py




