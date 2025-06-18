import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.pyplot as plt
from scipy.special import erf

import math
#<SOL:

# Definir las variables para la posición y velocidad del Sol en el sistema de coordenadas galácticas
x0 = 8.3 # kpc
y0 = 0 # kpc
z0 = 0.027 # kpc
vx0 = -11.1 # km/s
vy0 = 232.24 # km/s
vz0 = 7.25 # km/s


import numpy as np

#<SOL:

# Definir las variables para la posición y velocidad del Sol en el sistema de coordenadas galácticas
x0 = 8.3 # kpc
y0 = 0 # kpc
z0 = 0.027 # kpc
vx0 = -11.1 # km/s
vy0 = 232.24 # km/s
vz0 = 7.25 # km/s


# Calcular la velocidad tangencial del Sol en el plano galáctico
v0 = (vx0**2 + vy0**2)**0.5
d0 = math.sqrt(x0**2 + z0**2)
print(f"El Sol esta a {d0:.2f} kpc del centro de la galasia\n\
    y viaja a {v0:.2f} Km/seg tangencial mente")
#/SOL>

# Definimos los argumentos para el Sol
masa_entorno = 8.615e10  # Masa total de la Vía Láctea en unidades solares
masa_sol = 1.989e30  # Masa del Sol en kg
distancia_sol = 8  # Distancia del Sol al centro de la galaxia en kpc
v_max = 225  # Velocidad máxima de la galaxia en km/s


variables = 'source_id, ra, dec, parallax, radial_velocity, designation, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, pmra, pmdec'
n_ = ['source_id', 'ra', 'dec', 'parallax', 'radial_velocity', 'designation', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'pmra', 'pmdec']

def get_jobs_results(jobs):
    results = []
    for job in jobs:
        r = job.get_results()
        if r is not None:
            results.append(r)
    if len(results) > 0:
        return pd.concat(results, ignore_index=True)
    else:
        return None

def sql_gaia(n_estrellas):
    import os.path

    if os.path.isfile("datos_gaia.csv"):
        print("datos_gaia.csv existe no necesaria consulta a Gaia")
        return

    print("Haremos una consulta a Gaia para obtener datos")

    # Crear la consulta SQL
    sql = f"""
    SELECT TOP {n_estrellas} {variables}
    FROM gaiadr3.gaia_source
    """

    # Ejecutar la consulta y obtener los resultados
    print (sql)
    
    job = Gaia.launch_job_async(sql)
    result = job.get_results()
    # Convertimos los resultados a una tabla de astropy
    table = Table(result)

    # Escribimos la tabla en un archivo CSV
    table.write('datos_gaia.csv', format='csv', overwrite=True)

    return len(table)


def modelo_VN(masa_entorno, masa, radio, v_max=225):
    """
    Calcula la velocidad tangencial de una estrella utilizando el modelo newtoniano.

    Args:
    - masa_entorno (float): masa del entorno de la estrella en Masa Solar.
    - masa (float): masa de la estrella en Masa Solar.
    - radio (float): distancia de la estrella al centro de la galaxia en kpc.
    - v_max (float): velocidad de rotación máxima de la galaxia en km/s.

    Returns:
    - v_n (float): velocidad tangencial de la estrella en km/s.
    """
    G = 4.302e-6 # Constante gravitacional en unidades de kpc/Masa Solar (km/s)^2
    v_n = np.sqrt(G * masa_entorno / radio) * np.sqrt(1 + masa / masa_entorno) * (v_max / np.sqrt(radio))
    return v_n



def modelo_V1(masa_entorno, masa_estrella, radio, v_max=225):
    """
    Modelo Newtoniano para calcular la velocidad tangencial en el plano de la galaxia.
    
    Argumentos:
    - masa_entorno: masa total del entorno en unidades solares.
    - masa_estrella: masa de la estrella en unidades solares.
    - radio: radio de la órbita en kpc.
    - v_max: velocidad máxima en km/s.
    
    Retorna:
    - Velocidad tangencial en km/s.
    """
    G = 4.3e-6  # Constante gravitacional en unidades kpc^3/(M_sol*Gyr^2)
    masa_total = masa_entorno + masa_estrella
    v_c = np.sqrt(G * masa_total / radio)  # Velocidad circular en km/s
    v_n = v_max * np.sqrt(masa_entorno / masa_total)  # Velocidad de Navarro-Frenk-White en km/s
    v_t = np.sqrt(v_c**2 + v_n**2)  # Velocidad tangencial en km/s
    
    return v_t


def modelo_V2(masa_entorno, masa_estrella, radio, v_max=225):
    """
    Modelo de Einstein para calcular la velocidad tangencial en el plano de la galaxia.
    Incluye la propagación no instantánea de la fuerza gravitacional.
    
    Argumentos:
    - masa_entorno: masa total del entorno en unidades solares.
    - masa_estrella: masa de la estrella en unidades solares.
    - radio: radio de la órbita en kpc.
    - v_max: velocidad máxima de la galaxia (200-250) en km/s.
    Parametros:
    - c: velocidad de la luz en km/s.
    - beta: parámetro de ajuste para la propagación no instantánea de la fuerza gravitacional.
    
    Retorna:
    - Velocidad tangencial en km/s.
    """
    G = 4.3e-6  # Constante gravitacional en unidades kpc^3/(M_sol*Gyr^2)
    masa_total = masa_entorno + masa_estrella
    c = 299792458/1000  # C en Km/seg
    v_c = np.sqrt(G * masa_total / radio)  # Velocidad circular en km/s
    v_n = v_max * np.sqrt(masa_entorno / masa_total)  # Velocidad de Navarro-Frenk-White en km/s
    beta = 0.5  # Parámetro de ajuste para la propagación no instantánea de la fuerza gravitacional
    v_t = np.sqrt(v_c**2 + ((c**2 * v_n**2) / (v_n**2 + (beta * c**2))) * (1 - np.exp(-beta * radio)))  # Velocidad tangencial en km/s
    
    return v_t

# Atribución del matiz sobre la propagación no instantánea de la fuerza gravitacional a Víctor Estrada Díaz
# a partir de la sugerencia de Victor Estrada Diaz


sql_gaia(10000000)

# Crear dataframe vacío con columnas originales y nuevas 
#id,ra,dec,parallax,x,y,z,radio,distancia,vt,masa,vv
df = pd.DataFrame(columns=[\
    'id', 'ra', 'dec', 'paralax', \
    'x', 'y', 'z', 'radio', \
    'distancia', 'vt', 'masa', 'vv', \
    'masa_entorno', \
    'b', 'v_n', 'v_1', 'v_2'])
# Revisar el dataframe creado
print(df.head())

# unidades de distancia en Kpc
# unidades de masa en M(0) masas solares
# unidades de velocidad en Km/seg
# Identificacion y distancias respecto al Sol: 
#     'id', 'ra', 'dec', 'paralax'
# distancias respecto al centro galactico:
#     'x', 'y', 'z', 'radio'
# distancia al sol
#     'distancia'
# velocidad tangencial determinada de la observación
#     'vt', 
# velocidad radial dada por GAIA
#     'vv'
# ESTOS SON VALORES SIGUEN SE CALCuLARAN AHORA
# SON VALORES EXCLUIVAMENTE TEORICOS
# PARA CONTRASTAR CON LA VELOCIDAD DEDUCIDA 
# DE LOS DATOS DE GAIA (vt)
# Masa en el entorno de la estrella
#   'masa_entorno'
# Prametro (b) para las ecuaciones de campo de Einstein
#   'b'
# MODELO VN velocidad tangencial según Newton
#   'v_n'
# MODELO V1 velocidad tangencial según Einstein
# Este modelo no tiene en cuenta la velocidad finita
# de la propagación de la fuerza, toda deformación
# del espacio se propaga a velocidad C
# en el caso de las distancias en la galaxia esta 
# propagación no es algo a trivializar, este modelo V1
# no contempla esta velocidad C de propagación
#   'v_1'
# MODELO V2 velocidad tangencial corregida de Einstein
# donde si se tiene en cuenta la velocidad de propagacion
# de la deformacion del espacio
#   'v_2'

# Crear un DataFrame vacío con las columnas necesarias
df = pd.read_csv('estrellas.csv', usecols=range(12))
df['masa_entorno'] = 0
df['b'] = 0
df['v_n'] = 0
df['v_1'] = 0
df['v_2'] = 0

    
# Aqui ya leimos el csv y tenemos inicializado el trabajo

print(df.size)

masa_maxima = 1000
# Calcular la masa total de la muestra
df.loc[df['masa'] > masa_maxima, 'masa'] = np.nan
df = df.dropna(how="any", subset=["masa"])

masa_total = df['masa'].sum() * 10 #suma de la muestra
masa_via_lactea = 1.5 * 10**12  # masa de la Vía Láctea en masas solares
print("La masa de la Vía Láctea es aproximadamente:", masa_via_lactea, "masas solares")
#masa_total=masa_via_lactea

print(f"masa total:{masa_total:.0f} soles, muestras iniciales:{len(df) + len(df[df['masa']>masa_maxima])} validas:{len(df)}")

# Paso 1: Calcular la masa del entorno para cada estrella
df['masa_entorno'] = masa_total - df['masa']

# Aqui ya leimos el csv y tenemos inicializado el trabajo

print (df.size)



from scipy.stats import multivariate_normal
import numpy as np

# Obtener las medias y las desviaciones estándar de las distribuciones de masa y distancia
mu_masa, sigma_masa = df['masa'].mean(), df['masa'].std()
mu_distancia, sigma_distancia = df['radio'].mean(), df['radio'].std()

# Definir la función de densidad de probabilidad normal multidimensional
def multivariate_distribucion(x, y):
    rv = multivariate_normal([mu_masa, mu_distancia], 
                             [[sigma_masa**2, 0], [0, sigma_distancia**2]])
    return rv.pdf([x, y])

# Aplicar la función a cada par de valores de masa y distancia para obtener el array b
df['b'] = np.vectorize(multivariate_distribucion)(df['masa'], df['radio'])


i=0
for index, row in df.iterrows():
    
    masa=row.loc['masa']
    radio=row.loc['radio']
    masa_entorno=row.loc['masa_entorno']
    
    v_n = modelo_VN(masa_total,masa,radio,masa_entorno)
    df.at[index, 'v_n'] = v_n

    v_1 = modelo_V1(masa_entorno,masa,radio)
    df.at[index, 'v_1'] = v_1

    v_2 = modelo_V2(masa_entorno,masa,radio)
    df.at[index, 'v_2'] = v_2

#    i= i+1
#    if (i<10):
#        print (f"{index}  [{row.loc['b']:1f},{v_n:1f},{v_1:1f},{v_2:1f}]")
    

# Guardar el DataFrame resultante en un archivo csv
df.to_csv('estrellas_result.csv',header=True)



########REPRESENTACION########
# limitar df0 a las estrellas dentro del halo de la galaxia
# Definir el rango_halo como un porcentaje de max_radio
max_radio = 31
rango_halo = max_radio 
nucleo=0.5
print(f"Muestras iniciales {df.size}")
df = df.loc[df['radio']< rango_halo]
df = df.loc[df['radio']> nucleo]
print(f"Muestras seleccionadas {df.size}")
dd=df
dd_total = dd.size
dd = dd.loc [dd['radio']> 1.0]


import matplotlib.pyplot as plt
import os

os.remove('figura1.png')
os.remove('figura2.png')
os.remove('figura3.png')
os.remove('figura4.png')

# Gráfica 1
fig1, ax1 = plt.subplots(figsize=(6, 4))
ax1.scatter(df['radio'], df['vt'], s=1, c='green')
ax1.set_xlabel('Radio (pc)')
ax1.set_ylabel('FIG 1: Velocidad TANGENCIAL MEDIDA GAIA (km/s)')
ax1.scatter(d0, v0, c='red', marker='o', s=20)
fig1.savefig('figura1.png')
plt.close(fig1)

# Gráfica 2
fig2, ax2 = plt.subplots(figsize=(6, 4))
ax2.scatter(dd['radio'], dd['v_n'], s=1, c='orange')
ax2.set_ylabel('FIG 2: Velocidad TANGENCIAL SEGUN NEWTON (km/s)')
ax2.scatter(d0, v0, c='red', marker='o', s=20)
fig2.savefig('figura2.png')
plt.close(fig2)

# Gráfica 3
fig3, ax3 = plt.subplots(figsize=(6, 4))
ax3.scatter(dd['radio'], dd['v_1'], s=1, c='blue')
ax3.scatter(dd['radio'], dd['v_2'], s=1, c='green')
ax3.set_ylabel('FIG 3: Velocidad (azul) Modelo Einstein, (verde) Modelo VED(km/s)')
ax3.scatter(d0, v0, c='red', marker='o', s=20)
fig3.savefig('figura3.png')
plt.close(fig3)

# Gráfica 4
fig4, ax4 = plt.subplots(figsize=(6, 4))
ax4.scatter(dd['radio'], dd['b'], s=1, c='green')
ax4.set_ylabel('FIG 4: Distribución de la muestra')
ax4.scatter(d0, 0, c='red', marker='o', s=20)
fig4.savefig('figura4.png')
plt.show()
#plt.close(fig4)

from PIL import Image
import matplotlib.pyplot as plt

# Cargar las imágenes individuales
img1 = Image.open('figura1.png')
img2 = Image.open('figura2.png')
img3 = Image.open('figura3.png')
img4 = Image.open('figura4.png')

# Crear un cuadro para combinar las imágenes
cuadro = Image.new('RGB', (800, 600))

# Añadir las imágenes al cuadro
# Obtener el tamaño de cada imagen
width, height = img1.size

# Crear un cuadro grande para combinar las imágenes
cuadro = Image.new('RGB', (2 * width, 2 * height))

# Obtener el tamaño de cada imagen
width, height = img1.size
# Añadir las imágenes al cuadro
cuadro.paste(img1, (0, 0))
cuadro.paste(img2, (width, 0))
cuadro.paste(img3, (0, height))
cuadro.paste(img4, (width, height))

# Guardar el cuadro como una imagen
cuadro.save('cuadro.png')

# Mostrar el cuadro
plt.imshow(cuadro)
plt.axis('off')
plt.show()

print('FIN')