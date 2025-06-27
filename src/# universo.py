
from astroquery.sdss import SDSS
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.io import fits
from materia_oscura.src._ved_ import dis

"""Descripcion de la tabla en:
https://skyserver.sdss.org/dr12/en/help/browser/browser.aspx#&&history=description+Photoz+U
nombre	tipo	        largura	unidad	Ucd	        descripción	columnID
objID	bigint	        8	 	ID_MAIN	ID          único que apunta a la tabla GalaxyTag	1
z	    real	        4	 	REDSHIFT_PHOT	    corrimiento al rojo fotométrico; estimado por ajuste robusto a los vecinos más cercanos en un conjunto de referencia	2
zErr	real	        4	 	REDSHIFT_PHOT ERROR	error estimado del corrimiento al rojo fotométrico; si zErr=-9999, el ajuste falló, pero el corrimiento al rojo promedio de los vecinos aún podría estar disponible	3
nnCount	sint	        2	 	NÚMERO	            número de vecinos más cercanos después de excluir los valores atípicos; El valor máximo es 100, un valor mucho más pequeño indica una estimación deficiente	4
nnVol	real	        4	 	NÚMERO	            proporciona el volumen límite del espacio de color de los vecinos más cercanos de nnCount; Un valor grande indica una estimación deficiente	5
photoErrorClass	sint	2	 	CODE_MISC	        la clase de error fotométrico de la galaxia de 1 a 7, siendo 1 la mejor; un valor negativo muestra que el objeto a estimar está fuera de la caja del vecino más cercano k; zErr es una estimación confiable sólo si photoErrorClass=1	6
nnObjID	bigint	        8	 	ID_IDENTIFIER	    objID del (primer) vecino más cercano	7
nnSpecz	real	        4	 	CORRIMIENTO AL ROJO	corrimiento al rojo espectroscópico del (primer) vecino más cercano	8
nnFarObjID	bigint	    8	 	ID_IDENTIFIER	    objID del vecino más lejano	9
nnAvgZ	real	        4	 	CORRIMIENTO AL ROJO	corrimiento al rojo promedio de los vecinos más cercanos; Si es significativamente diferente de Z, esta podría ser una mejor estimación que Z.	10
distMod	real	        4	 	PHOT_DIST-MOD	el módulo de distancia para Omega=0.2739, Lambda=0.726, h=0.705 cosmología	11
lumDist	real	        4	 	PHOT_LUM-DIST	la distancia de luminosidad en Mpc para Omega=0.2739, Lambda=0.726, h=0.705 cosmología	12
Chisq	real	        4	 	STAT_LIKELIHOOD	El valor chi-cuadrado para el ajuste mínimo de la plantilla chi-cuadrado (no reducido, 4 grados de libertad)	13
rnorm	real	        4	 	STAT_MISC	el valor de la norma euclidiana residual para el ajuste mínimo de la plantilla chi-cuadrado	14
bestFitTemplateID	Int	4	 	ID_IDENTIFIER	identificador de la plantilla más adecuada; si bestFitTemplateID=0, todas las columnas siguientes no son válidas y pueden rellenarse con el valor -9999	15
sintetizador	real	4	 	PHOT_SYNTH-MAG PHOT_SDSS_U	Magnitud U sintética calculada a partir de la plantilla ajustada	16
synthG	real	        4	 	PHOT_SYNTH-MAG PHOT_SDSS_G	Magnitud g sintética calculada a partir de la plantilla ajustada	17
sintetizador	real	4	 	PHOT_SYNTH-MAG PHOT_SDSS_R	Magnitud R' sintética calculada a partir de la plantilla ajustada	18
synthI	real	        4	 	PHOT_SYNTH-MAG PHOT_SDSS_I	Magnitud I' sintética calculada a partir de la plantilla ajustada	19
sintetizador	real	4	 	PHOT_SYNTH-MAG PHOT_SDSS_Z	Magnitud z sintética calculada a partir de la plantilla ajustada	20
kcorrU	real	        4	 	PHOT_K-CORRECCIÓN	k corrección para z=0	21
kcorrG	real	        4	 	PHOT_K-CORRECCIÓN	k corrección para z=0	22
kcorrR	real	        4	 	PHOT_K-CORRECCIÓN	k corrección para z=0	23
kcorrI	real	        4	 	PHOT_K-CORRECCIÓN	k corrección para z=0	24
kcorrZ	real	        4	 	PHOT_K-CORRECCIÓN	k corrección para z=0	25
kcorrU01	real	    4	 	PHOT_K-CORRECCIÓN	k corrección para z=0,1	26
kcorrG01	real	    4	 	PHOT_K-CORRECCIÓN	k corrección para z=0,1	27
kcorrR01	real	    4	 	PHOT_K-CORRECCIÓN	k corrección para z=0,1	28
kcorrI01	real	    4	 	PHOT_K-CORRECCIÓN	k corrección para z=0,1	29
kcorrZ01	real	    4	 	PHOT_K-CORRECCIÓN	k corrección para z=0,1	30
absMagU	real	        4	 	PHOT_ABS-MAG PHOT_SDSS_U	Resto de la magnitud de los abdominales	31
absMagG	real	        4	 	PHOT_ABS-MAG PHOT_SDSS_G	Resto Marco G' ABS magnitud	32
absMagR	real	        4	 	PHOT_ABS-MAG PHOT_SDSS_R	Resto Frame R' ABS magnitud	33
absMagI	real	        4	 	PHOT_ABS-MAG PHOT_SDSS_I	Resto Marco I' Magnitud ABS	34
absMagZ	real	        4	 	PHOT_ABS-MAG PHOT_SDSS_Z	Resto del marco Z' Magnitud del ABS	35
"""
# Número de objetos que deseas obtener
num_objects = 10000

# Consulta al SDSS para obtener una muestra aleatoria de objetos por distancia
sql1 = f"""
    SELECT TOP {num_objects}
         objID, z, class  
    FROM
        SpecOb

    """
sql2 = f"""
    SELECT TOP {num_objects}
        objID,z,distMod,lumDist    
    FROM
        2dFGRS
    WHERE
        z <= 30.0
    """


sql3 = f"""SELECT TOP {num_objects} 
   objID, z
FROM 2dFGRS
WHERE z <= 1
""" 
query = sql2

try:
    result = SDSS.query_sql(query)
except Exception as ex:
    print(f"Error: {ex}")


data = []  # Lista para almacenar los datos
for row in result:
    data.append([row['objID'],row['z'], row['distMod'], row['lumDist']])
    
    
columns = ['objID','z', 'distMod', 'lumDist']  # Nombres de las columnas
f = pd.DataFrame(data, columns=columns)  # Crear el DataFrame
# Filtrar los valores menores o iguales a -999.0 en z y distMod
df = f[(f['z'] > -999.0) & (f['distMod'] > -999.0)]

print(f"maximos dismod {df['distMod'].max()} lumDist {df['lumDist'].max()} z {df['z'].max()}")
print(f"minimos dismod {df['distMod'].min()} lumDist {df['lumDist'].min()} z {df['z'].min()}")

ng = 1
figura = f'galaxias{ng}'
dis(
    figura,
    df['distMod'].values, # Módulo de distancia (distMod)
    df['z'].values,       # Corrimiento al rojo (z)
    'Módulo de distancia vs. Corrimiento al rojo', # título
    'DistMod',            # subtitulo (puedes dejarlo en blanco si no lo necesitas)
    'blue',               # color de los puntos
    'Módulo de distancia',      # etiqueta del eje y
    'Corrimiento al rojo (z)',  # etiqueta del eje x
    0, 0                  # xlim y ylim (no se aplican límites en este caso)
)
ng += 1
figura = f'galaxias{ng}'
dis(
    figura,
    df['lumDist'].values,   # Distancia por luminosidad (Lumdist)       
    df['z'].values,         # Corrimiento al rojo (z)
    'Gráfico de Corrimiento al rojo vs. Módulo de distancia', # título
    'Por Luminosidad',      # subtitulo (puedes dejarlo en blanco si no lo necesitas)
    'blue',                 # color de los puntos
    'distancia por luminosidad', # etiqueta del eje x
    'Corrimiento al rojo (z)',   # etiqueta del eje y
    0, 0                         # xlim y ylim (no se aplican límites en este caso)
)

# Convertir corrimiento al rojo (z) en tiempo en millones de años (conversión aproximada)
# Supongamos que se utiliza una constante de Hubble de 70 km/s/Mpc
hubble_constant = 70  # km/s/Mpc
light_speed = 299792.458  # km/s (velocidad de la luz)
mega_years_to_years = 1000000  # conversión de mega años a años

def redshift_to_time(z):
    return (-1 / hubble_constant) * (light_speed / mega_years_to_years) * z

# Aplicar la función de conversión a la columna 'z' para obtener el tiempo en millones de años
df['tiempo'] = df['z'].apply(redshift_to_time)

# Crear la gráfica
ng += 1
figura = f'galaxias{ng}'
dis(
    figura,
    df['tiempo'].values,  # Usar el tiempo en el eje x
    df['distMod'].values,
    'Gráfico de Tiempo vs. Módulo de distancia',
    'Tiempo - DistMod',
    'blue',
    'Tiempo (millones de años desde el Big Bang)',
    'Módulo de distancia',
    0, 0
)
# ... (código previo)

# Convertir tiempo en millones de años a corrimiento al rojo (z)
def time_to_redshift(tiempo):
    return (-hubble_constant * tiempo * mega_years_to_years) / light_speed

# Crear un rango de tiempo desde -14700 millones de años hasta 0 (presente)
tiempo_range = np.linspace(-14700, 0, num=1000)

# Calcular el corrimiento al rojo correspondiente para cada valor de tiempo en el rango
redshift_values = [time_to_redshift(tiempo) for tiempo in tiempo_range]

# Crear la gráfica de corrimiento al rojo vs. tiempo
ng += 1
figura = f'galaxias{ng}'
dis(
    figura,
    tiempo_range,         # Usar el tiempo en el eje x
    redshift_values,      # Usar los valores de corrimiento al rojo en el eje y
    'Evolución del Corrimiento al Rojo desde el Big Bang',
    'Tiempo - Corrimiento al rojo',
    'blue',
    'Tiempo (millones de años desde el Big Bang)',
    'Corrimiento al Rojo (z)',
    0, 0
)

# ...
import requests
from astropy.io import fits

# URL desde la cual se descargará el archivo FITS
fits_url = 'URL_DEL_ARCHIVO_FITS_AQUI'

# Ruta donde se guardará el archivo FITS descargado
downloaded_filename = 'spAll-v5_10_0.fits'

# Descargar el archivo FITS
response = requests.get(fits_url)
with open(downloaded_filename, 'wb') as f:
    f.write(response.content)

# Abrir el archivo FITS descargado
with fits.open(downloaded_filename) as fits_file:
    # Obtiene el encabezado de la primera extensión
    header = fits_file[1].header

    # Imprime los nombres de todas las claves disponibles en el encabezado
    print(header.keys())

    # Lee los datos de las columnas especificadas en 'columns'
    columns = ['PLATE', 'MJD', 'FIBERID', 'Z', 'ZWARNING', 'Z_ERR']
    data = fits_file[1].data

    # Crea un diccionario para almacenar las columnas seleccionadas
    selected_data = {col: data[col] for col in columns}

# 'selected_data' contiene las columnas seleccionadas
print(selected_data)