import pandas as pd
import numpy as np
from astroquery.gaia import Gaia
from astropy.table import Table
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from scipy.interpolate import make_interp_spline
import  math
import os.path

variables = 'source_id, ra, dec, parallax, radial_velocity, designation, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, pmra, pmdec'
n_ = variables.split(', ')

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

def calcular_teff(g, bp, rp):
    """
    Esta función calcula la temperatura efectiva (teff) 
    utilizando las magnitudes en banda G, BP y RP
    """
    # Calcular el color (BP-RP)
    bp_rp = bp - rp
    
    # Coeficientes de las relaciones empíricas
    c1 = 3.939
    c2 = 0.6536
    c3 = 0.02305
    c4 = 0.004012
    
    # Calcular la temperatura efectiva (teff)
    teff = c1 + c2 * bp_rp + c3 * bp_rp**2 + c4 * bp_rp**3
    teff = 10**(teff)
    
    # Corregir la temperatura efectiva por la extinción
    ag = 0.859 # Extinción en banda G
    teff = teff / (1 - ag/3.1)
    
    return teff

import math

def calcular_coordenadas_galacticas(ra, dec, parallax):
    # Distancia en parsecs desde el Sol a la estrella
    distancia_parsecs = 1 / (float(parallax) * 0.001)
    
    # Coordenadas cartesianas de la estrella respecto al Sol en parsecs
    x_parsecs = distancia_parsecs * math.cos(float(dec)) * math.cos(float(ra))
    y_parsecs = distancia_parsecs * math.cos(float(dec)) * math.sin(float(ra))
    z_parsecs = distancia_parsecs * math.sin(float(dec))
    
    # Coordenadas cartesianas de la estrella respecto al centro galáctico en kiloparsecs
    x_galacticas = (x_parsecs - 8.0) / 1000.0
    y_galacticas = y_parsecs / 1000.0
    z_galacticas = z_parsecs / 1000.0
    
    return x_galacticas, y_galacticas, z_galacticas


def calcular_masa(teff, Gmag, distancia, color):
    Mbol_sun = -26.832
    BC = 0.05
    try:
        M_G = Gmag - 5 * math.log10(distancia/10)
        Mbol = M_G + 0.006 * (teff - 5777)
        M = 10 ** ((Mbol - Mbol_sun - BC) / 2.5)
    except:
        M = np.exp(1.43 + 1.66 * color)
    return M

import math

def calcular_velocidad(pmra, pmdec, distancia):
    # Calcular el movimiento propio total
    movimiento_propio = math.sqrt(pmra**2 + pmdec**2)
    
    # Convertir el movimiento propio de arcsec/año a km/seg
    # La constante 4.74047 se utiliza para la conversión
    movimiento_propio_km_seg = movimiento_propio * distancia * 4.74047
    
    # Calcular la componente en z de la coordenada Galáctica (velocidad en la dirección vertical)
    velocidad_z = movimiento_propio_km_seg * math.sin(math.radians(90 - math.degrees(math.atan2(pmdec, pmra))))
    
    return movimiento_propio_km_seg, velocidad_z

# Ejemplo de uso
pmra = 2.5  # Componente de movimiento propio en ascensión recta en arcsec/año
pmdec = -1.8  # Componente de movimiento propio en declinación en arcsec/año
distancia = 1000  # Distancia en parsecs

movimiento_propio_km_seg, velocidad_z = calcular_velocidad(pmra, pmdec, distancia)

print("Movimiento propio: {:.2f} km/seg".format(movimiento_propio_km_seg))
print("Velocidad en la dirección z: {:.2f} km/seg".format(velocidad_z))


def estrellas():
    df = pd.read_csv("datos_gaia.csv")
    print (f"filas de consulta SQL: {df.size}")
    with open("estrellas.csv", "w") as f:
        n=0
        f.write("id,ra,dec,parallax,x,y,z,radio,distancia,vt,masa,vv\n")
        df_clean = df.dropna(subset=['pmra', 'pmdec']) # elimina filas con NaN en pmra o pmdec
        print (f"filtrado a datos validos: {df_clean.size}")

        for line in df_clean.values:
            n=n+1
            id = line[n_.index('source_id')]
            ra = line[n_.index('ra')]
            dec = line[n_.index('dec')]
            parallax = line[n_.index('parallax')]
            distancia = 1/float(parallax)
            x, y, z = calcular_coordenadas_galacticas(ra, dec, parallax)
            radio = math.sqrt(x**2 + y**2 + z**2)
            pmra = line[n_.index('pmra')]
            pmdec = line[n_.index('pmdec')]
            I=math.atan(x/y)
            B=math.atan(z/math.sqrt(x**2+y**2))
            vt = calcular_velocidad(pmra, pmdec, distancia)
            g = line[n_.index('phot_g_mean_mag')]
            bp = line[n_.index('phot_bp_mean_mag')]
            rp = line[n_.index('phot_rp_mean_mag')]
            teff = calcular_teff(g, bp, rp)
            masa = calcular_masa(teff, g, distancia,bp-rp)
            vv = line[n_.index('radial_velocity')]
            #if radio > 1.2 and radio < 40:
            f.write(f"{id},{ra},{dec},{parallax},{x},{y},{z},{radio},{distancia},{vt},{masa},{vv}\n")
            n=n+1
    print (f"Estrellas con datos validos:{n}")

sql_gaia(1000000)
estrellas()

# cargar los datos del archivo estrellas.csv
df0=pd.read_csv('estrellas.csv', engine='python')


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

# limitar df0 a las estrellas dentro del halo de la galaxia
# Definir el rango_halo como un porcentaje de max_radio
max_radio = 31
rango_halo = max_radio 
#obtenido del histograma 95% de las estrellas

nucleo=0.59
print(f"Muestras iniciales {df0.size}")
df = df0.loc[df0['radio']< rango_halo]
dd = df.loc[df['radio']> nucleo]
print(f"Muestras seleccionadas {df.size}")

dd.dropna(subset=['vv']) # elimina filas con NaN
dd_total = dd.size
dd = dd.loc [dd['radio']> 0]
dd = dd.loc [dd['radio']<rango_halo]
 
print(f"Muestras con velocidad radial valida{dd.size}")

n, bins = np.histogram(df['radio'], bins=10)

# calcular la distribución acumulada
cumulative_n = np.cumsum(n)

# encontrar la distancia al percentil 95
percentile_95_index = np.argmax(cumulative_n > 0.95 * len(df))
distancia_percentil_95 = (bins[percentile_95_index] + bins[percentile_95_index + 1]) / 2

import os
archivos = ["figura01.png", "figura02.png", "figura03.png","figura04.png"]

for archivo in archivos:
    try:
        os.remove(archivo)
    except FileNotFoundError:
        pass
    
# Gráfica 1
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='3d')

# Gráfico 3D de las coordenadas x, y, z
ax1.scatter(dd['x'], dd['y'], dd['z'], color='black', s=0.5)
ax1.set_xlabel('X (pc)')
ax1.set_ylabel('Y (pc)')
ax1.set_zlabel('Z (pc)')
#el sol en rojo
#ax1.scatter(x0, y0, z0, c='r', marker='o', s=2)

fig1.savefig('figura01.png')
plt.show()
plt.close(fig1)

# gráfico de radio y velocidad tangencial (vt)
# Gráfica 2
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.text(0.5, 1.05, f'Número de estrellas: {len(df)}', ha='center', va='center', transform=ax2.transAxes)
ax2.scatter(df['radio'], df['vt'])
ax2.set_xlabel('Radio (pc)')
ax2.set_ylabel('Velocidad tangencial (km/s)')
#el sol en rojo
ax2.scatter(d0, v0, c='r', marker='o', s=20)
fig2.savefig('figura02.png')
plt.show()
plt.close(fig2)

# Gráfica 3
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
# plotear el histograma y la línea del percentil 95
sns.kdeplot(df['radio'], ax=ax3, color='blue', linewidth=2, alpha=0.5)
ax3.axvline(distancia_percentil_95, color='red', linestyle='--', linewidth=2)
ax3.set_ylabel('Densidad de probabilidad')
ax3.set_title(f'Percentil 95 = {distancia_percentil_95:.2f}')
fig3.savefig('figura03.png')
plt.close(fig3)

# Gráfica 4
fig4, ax4 = plt.subplots(figsize=(6, 4))
# gráfico de radio y velocidad radial (vt)
ax4.scatter(dd['radio'], dd['vv'], s=4)
ax4.set_xlabel('Radio (pc)')
ax4.set_ylabel('Velocidad radial (km/s)')
ax4.set_title(f'Estrellas con v. radial = {dd.size} de {dd_total}')
ax4.scatter(d0, v0, c='r', marker='o', s=20)
fig4.savefig('figura04.png')
plt.close(fig4)

from PIL import Image
import matplotlib.pyplot as plt

# Cargar las imágenes individuales
img1 = Image.open('figura01.png')
img2 = Image.open('figura02.png')
img3 = Image.open('figura03.png')
img4 = Image.open('figura04.png')


print("FIN")