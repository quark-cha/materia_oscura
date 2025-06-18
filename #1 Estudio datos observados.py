import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import make_interp_spline
import  math
import numpy as np
import os.path
from scipy.stats import multivariate_normal
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord, Galactic, ICRS
import astropy.units as u
import pandas as pd

import networkx as nx

global recalcular 
recalcular = False

### Importante ###
# Borrar datos_gaia si desea recargar datos desde el proyecto GAIA
# Si cambia recalcular a:
#       True      se rehará ESTRELLAS_RESULT.CSV con nuevos clculos
#       False     se omitiran los calculos y se cargaran de ESTRELLAS_RESULT.CSV

### GLOBALES ###

global plt
global m_iniciales #muestras iniciales
global m_situacion #muestras situacio apropiada
global m_velocidad #muestras velocidad valida
global m_masa      #muestras con masa valida
global masa_total 
global x0,y0,z0,vx0,vy0,vz0,v0,d0

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
print(f"El Sol esta a {d0:.2f} kpc del centro de la galaxia\n\
  y viaja a {v0:.2f} Km/seg tangencial mente")

### QUERY SQL ###
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

### CONSULTA A GAIA ###
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

### CALCULOS PREVIOS ###

def vector_v(ra, dec, parallax, pm_ra, pm_dec, radial_velocity, test=False):
    # Obtén los valores observados de la estrella
    distance_pc = 1.0 / (parallax)  # Distancia de la estrella en parsecs
    if distance_pc < 0: distance_pc = -distance_pc

    # Factor de conversión de mas/yr a km/s para distancias en parsecs
    conversion_factor = 4.74

    # Convierte las componentes de velocidad propia de 'mas/yr' a 'km/s'
    pm_ra_km_s = pm_ra * distance_pc * conversion_factor
    pm_dec_km_s = pm_dec * distance_pc * conversion_factor

    # Convertir los ángulos ra y dec a radianes
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    
    # Calcular las componentes x, y, z del vector de movimiento propio en km/s
    vx = pm_ra_km_s * np.cos(ra_rad)
    vy = pm_ra_km_s * np.sin(ra_rad) + pm_dec_km_s * np.cos(ra_rad)
    vz = pm_dec_km_s * np.sin(ra_rad)
    
    # Estas velocidades son las que observa el sol y por tanto no tienen en cuenta su propio movimiento
    # Necesitamos restarle las velocidades del propio movimiento del sol
    
    # Velocidad del Sol en el marco GSR en km/s
    v_sun_x = 11.1  # Componente x de la velocidad del Sol en km/s
    v_sun_y = 232.24  # Componente y de la velocidad del Sol en km/s
    v_sun_z = 7.25  # Componente z de la velocidad del Sol en km/s

    # Restar vectorialmente las velocidades del propio movimiento del Sol
    vx_real = vx - v_sun_x
    vy_real = vy - v_sun_y
    vz_real = vz - v_sun_z
    
    # Agregar la velocidad radial si está disponible
    if radial_velocity is not None:
        vz_real += radial_velocity

    # Calcular módulos
    mod_v = np.sqrt(vx_real**2 + vy_real**2 + vz_real**2)
    mod_t = np.sqrt(vx_real**2 + vy_real**2)

    # Calcular divergencia
    divergencia = mod_v - mod_t

    # Crear DataFrame
    data = {
        'vx': [vx_real],
        'vy': [vy_real],
        'vz': [vz_real],
        'mod_v': [mod_v],
        'mod_t': [mod_t],
        'divergencia': [divergencia]
    }

    df = pd.DataFrame(data)
    return df

def rv_to_gsr(c, v_sun=None):
    """Transform a barycentric radial velocity to the Galactic Standard of Rest
    (GSR).

    The input radial velocity must be passed in as a

    Parameters
    ----------
    c : `~astropy.coordinates.BaseCoordinateFrame` subclass instance
        The radial velocity, associated with a sky coordinates, to be
        transformed.
    v_sun : `~astropy.units.Quantity`, optional
        The 3D velocity of the solar system barycenter in the GSR frame.
        Defaults to the same solar motion as in the
        `~astropy.coordinates.Galactocentric` frame.

    Returns
    -------
    v_gsr : `~astropy.units.Quantity`
        The input radial velocity transformed to a GSR frame.

    """
    if v_sun is None:
        v_sun = coord.Galactocentric().galcen_v_sun.to_cartesian()

    gal = c.transform_to(coord.Galactic)
    cart_data = gal.data.to_cartesian()
    unit_vector = cart_data / cart_data.norm()
    v_proj = v_sun.dot(unit_vector)
    return c.radial_velocity + v_proj


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

def estrellas():
    df = pd.read_csv("datos_gaia.csv")
    num_estrellas_sin_vr = df[df['radial_velocity'].isnull() | (df['radial_velocity'] == 0)].shape[0]
    print(f"Número de estrellas sin velocidad radial: {num_estrellas_sin_vr}")
    df['radial_velocity'] = df['radial_velocity'].fillna(0)
    print (f"filas de consulta SQL: {df.size}")
    n=0
    with open("estrellas.csv", "w") as f:
        f.write("id,ra,dec,parallax,x,y,z,radio,distancia,pmra,pmdec,vt,masa,vv\n")
        df_clean = df.dropna(subset=['pmra', 'pmdec']) # elimina filas con NaN en pmra o pmdec
        print (f"filtrado a datos validos: {df_clean.size}")
        for line in df_clean.values:
            n=n+1    
            id = line[n_.index('source_id')]
            ra = line[n_.index('ra')]
            dec = line[n_.index('dec')]
            parallax = line[n_.index('parallax')]
            distancia = 1/(float(parallax)*0.001)
            x = distancia * math.cos(float(ra)) * math.cos(float(dec)) /1000.0
            y = distancia * math.sin(float(ra)) * math.cos(float(dec)) /1000.0
            z = distancia * math.sin(float(dec)) /1000.0
            radio = math.sqrt(x**2 + y**2 + z**2)
            pmra = line[n_.index('pmra')]
            pmdec = line[n_.index('pmdec')]
            I=math.atan(x/y)
            B=math.atan(z/math.sqrt(x**2+y**2))
            vt = 4.74 * radio * math.sqrt(pmra**2 + (pmdec - 4.95*math.cos(I)*math.cos(B)- 0.169*math.sin(I)*math.cos(B) + 0.02*math.sin(B))**2)
            g = line[n_.index('phot_g_mean_mag')]
            bp = line[n_.index('phot_bp_mean_mag')]
            rp = line[n_.index('phot_rp_mean_mag')]
            teff = calcular_teff(g, bp, rp)
            masa = calcular_masa(teff, g, distancia,bp-rp)
            vv = line[n_.index('radial_velocity')]
            #if radio > 1.2 and radio < 40:
            f.write(f"{id},{ra},{dec},{parallax},{x},{y},{z},{radio},{distancia},{pmra},{pmdec},{vt},{masa},{vv}\n")

    print (f"Estrellas con datos validos:{n}")

def multivariate_distribucion(x, y):     
    rv = multivariate_normal([mu_masa, mu_distancia],\
        [[sigma_masa**2, 0], [0, sigma_distancia**2]])
    return rv.pdf([x, y])

def leer_estrellas():
    # Crear dataframe vacío con columnas originales y nuevas 
    #id,ra,dec,parallax,x,y,z,radio,distancia,vt,masa,vv,masa_entorno,b,v_n,v_1,v_2,v_e,v_3,v_4,v_5,v_6,v_7
    df = pd.DataFrame(columns=[\
        'id', 'ra', 'dec', 'paralax', \
        'x','y','z','radio', \
        'vt',\
        'masa','vv','masa_entorno', 'b', \
        'v_n', 'v_1', 'v_2', 'v_e', 'v_3', 'v_4', 'v_5', 'v_6', 'v_7'])

    # Leer el archivo CSV en un nuevo DataFrame
    df_csv = pd.read_csv('estrellas.csv')

    # Combinar los DataFrames asegurándose de que las columnas coincidan
    df = pd.concat([df_csv,df], sort=False)

    # Inicializar las columnas no utilizadas con ceros
    columnas_faltantes = set(df.columns).difference(df_csv.columns)
    for columna in columnas_faltantes:
        df[columna] = 0
    
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
    # velocidad tangencial determinada de la observación no GSM
    #     'vt', 
    # velocidad radial dada por GAIA
    #     'vv'
    # ESTOS SON VALORES SIGUEN SE CALCULARAN AHORA
    # SON VALORES EXCLUIVAMENTE TEORICOS
    # PARA CONTRASTAR CON LA VELOCIDAD DEDUCIDA 
    # DE LOS DATOS DE GAIA (vt)
    
    print("filas leidas ... ",df.size)

    ### CALCULO de distribucion de MASAS df:['b'] ###
    
    masa_maxima = 1000

    # Calcular la masa total de la muestra

    df.loc[df['masa'] > masa_maxima, 'masa'] = np.nan
    df = df.dropna(how="any", subset=["masa"])

    global masa_total 
    masa_total = df['masa'].sum() * 10 #suma de la muestra
    masa_via_lactea = 1.5 * 10**12  # masa de la Vía Láctea en masas solares
    print("La masa de la Vía Láctea es aproximadamente:", masa_via_lactea, "masas solares")
    
    #masa_total=masa_via_lactea
    print(f"masa total:{masa_total:.0f} soles, muestras iniciales:{len(df) + len(df[df['masa']>masa_maxima])} validas:{len(df)}")

    # Paso 1: Calcular la masa del entorno para cada estrella
    df['masa_entorno'] = masa_total - df['masa']
    print (df.size)
    return df


def completar_estrellas(df,modelos=True):

    # Obtener las medias y las desviaciones estándar de las distribuciones de masa y distancia
    global mu_masa, mu_distancia, sigma_masa, sigma_distancia
    mu_masa, sigma_masa = df['masa'].mean(), df['masa'].std()
    mu_distancia, sigma_distancia = df['radio'].mean(), df['radio'].std()

    # Definir la función de densidad de probabilidad normal multidimensional

    # Aplicar la función a cada par de valores de masa y distancia para obtener el array b
    df['b'] = np.vectorize(multivariate_distribucion)(df['masa'], df['radio'])

    PRE_Modelo_VED(df)
    df['v_7'] = velocidad_propia

    df['v_g'] = PRE_GRAFO(df)
    
    
    if modelos: 
        for index, row in df.iterrows():
            masa=row.loc['masa']
            radio=row.loc['radio']
            masa_entorno=row.loc['masa_entorno']
            probabilidad=row.loc['b']
            radial_velocity=row.loc['vv']
            ra=row.loc['ra']
            dec=row.loc['dec']
            parallax=row.loc['parallax']
            pmra=row.loc['pmra']
            pmdec=row.loc['pmdec']
        
            ### Obtiene los datos de los modelos para df

            v_n = modelo_VN(masa_total,masa,radio,masa_entorno)
            df.at[index, 'v_n'] = v_n
            v_1 = modelo_V1(masa_entorno,masa,radio)
            df.at[index, 'v_1'] = v_1
            v_2 = modelo_V2(masa_entorno,masa,radio)
            df.at[index, 'v_2'] = v_2
            v_3 = modelo_V3(masa_total, masa, radio)
            df.at[index, 'v_3'] = v_3
            v_4 = modelo_V4(masa_total,masa,radio)
            df.at[index, 'v_4'] = v_4
        
            # Calcular la velocidad tangencial desde la Galaxia
            v=vector_v(ra, dec, parallax, pmra, pmdec, radial_velocity)
            #'vx': 'vy': 'vz': 'mod_v': 'mod_t': 'divergencia'
            v_5 = v['mod_v'].values[0]
            v_6 = v['mod_t'].values[0]
        
            df.at[index, 'v_5'] = v_5
            df.at[index, 'v_6'] = v_6
        
        
        
    ### Fin de asignacion de datos de modelos    


    # limitar df a las estrellas dentro del halo de la galaxia

    # Definir el rango_halo como un porcentaje de max_radio

    max_radio = 40

    rango_halo = max_radio 

    #obtenido del histograma 95% de las estrellas


    nucleo=0.01

    print(f"Muestras iniciales {df.size}")

    m_iniciales = df.size

    dd = df.loc[df['radio']< rango_halo]

    dd = dd.loc[dd['radio']> nucleo]

    print(f"Por situacion seleccionadas {dd.size}")

    m_situacion= dd.size

    # elimina filas con NaN

    dd.dropna(subset=['vv'])

    dd.dropna(subset=['masa'])

    dd=dd.loc[dd['masa']<1200]

    dd=dd.loc[dd['masa']>0.1]

    print(f"Con masa valida seleccionadas 2 {dd.size}")

    m_valida= dd.size

    dd = dd.loc [dd['radio']> 0]

    dd = dd.loc [dd['radio']<rango_halo]

    m_velocidad = dd.size

    print(f"Con velocidad radial valida {dd.size}")


    n, bins = np.histogram(df['radio'], bins=10)


    # calcular la distribución acumulada
    cumulative_n = np.cumsum(n)


    # encontrar la distancia al percentil 95

    percentile_95_index = np.argmax(cumulative_n > 0.95 * len(df))

    distancia_percentil_95 = (bins[percentile_95_index] + bins[percentile_95_index + 1]) / 2

    return dd




######### MODELOS #########

# Cálculo vectorial del premodelo VED
import math
import numpy as np
import networkx as nx

def interaccion_gravitatoria_significativa(x1, z1, masa1, x2, z2, masa2, umbral):
    distancia = math.sqrt((abs(x1 - x2))**2 + (abs(z1 - z2))**2)
    G = 6.67430e-11  # Constante gravitacional
    fuerza = G * (masa1 * masa2) / distancia ** 2
    return fuerza < umbral

def PRE_GRAFO(df):
    umbral = 1e-9  # Valor de umbral para la fuerza gravitatoria (ajústalo según tus necesidades)
    nodos = np.arange(len(df))
    G = nx.Graph()
    G.add_nodes_from(nodos)
    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            row1 = df.loc[i]
            row2 = df.loc[j]
            if interaccion_gravitatoria_significativa(row1['x'], row1['z'], row1['masa'],
                                                     row2['x'], row2['z'], row2['masa'], umbral):
                G.add_edge(i, j)

    shortest_paths = nx.shortest_path(G)
    communities = nx.algorithms.community.greedy_modularity_communities(G)
    degree_centrality = nx.degree_centrality(G)

    velocidad_resultante = {}
    for nodo in G.nodes:
        vecinos = list(G.neighbors(nodo))
        masa_nodo = df.loc[nodo, 'masa']
        suma_velocidades = df.loc[vecinos, 'velocidad'].sum()
        velocidad_resultante[nodo] = suma_velocidades / masa_nodo

    return velocidad_resultante


def PRE_Modelo_VED(df):
    global velocidad_propia
        
    # Convertir masa de estrella a kilogramos
    M_estrella = df['masa'] * 1.989e30  # 1 masa solar = 1.989 × 10^30 kg

    # Convertir radio a metros
    radio_m = df['radio'] * 3.086e16  # 1 kiloparsec = 3.086 × 10^16 metros

    d_b = df['b']
    print(f"B mínima: {df['b'].min()}, máxima: {df['b'].max()}")

    c = 299792458.0  # Velocidad de la luz en metros por segundo

    # Calcular la fuerza resultante del efecto propagativo para todas las estrellas en newtons
    F_P = c**2 * M_estrella * radio_m * d_b
    print(f"F_P mínima: {F_P.min()}, máxima: {F_P.max()}")

    G = 6.67430e-11  # Constante gravitacional en m^3 kg^-1 s^-2

    # Calcular la fuerza resultante de la distribución de masa para todas las estrellas en newtons
    F_D = G * M_estrella / radio_m ** 2
    print(f"F_D mínima: {F_D.min()}, máxima: {F_D.max()}")

    # Calcular la fuerza resultante total para todas las estrellas en newtons
    F_resultante = F_D - F_P

    # Calcular la velocidad tangencial
    velocidad_tangencial = np.sqrt(np.maximum(F_resultante * radio_m / M_estrella, 0))
    #velocidad_tangencial = np.sqrt(F_resultante * radio_m / M_estrella)

    # Convertir la velocidad
    # tangencial a kilómetros por segundo
    velocidad_tangencial_kmps = velocidad_tangencial / 1000

    print(f"Velocidad tangencial mínima: {velocidad_tangencial_kmps.min()}, Velocidad tangencial máxima: {velocidad_tangencial_kmps.max()}")
 
    velocidad_propia = velocidad_tangencial_kmps

### VN: df['V_N']
def modelo_VN(masa_entorno=None, masa=None, radio=None, v_max=225, describe=None):
    """
    MODELO VN: Modelo de Newton
    Calcula la velocidad tangencial de una estrella utilizando el modelo newtoniano.
    Args:
    - masa_entorno (float): masa del entorno de la estrella en masas solares.
    - masa (float): masa de la estrella en masas solares.
    - radio (float): distancia de la estrella al centro de la galaxia en kiloparsecs.
    - v_max (float): velocidad de rotación máxima de la galaxia en km/s.
    Returns:
    - v_n (float): velocidad tangencial de la estrella en km/s.
    """
    if describe:
        print(modelo_VN.__doc__)
        return

    # Conversión de unidades
    masa_entorno_kg = masa_entorno * 1.988e30  # Conversión de masas solares a kilogramos
    masa_kg = masa * 1.988e30
    radio_km = radio * 3.086e13  # Conversión de kiloparsecs a kilómetros
    G = 6.67430e-20  # Constante gravitacional en unidades de km^3/(kg*s^2)
    v_n = np.sqrt(G * masa_entorno_kg / radio_km) * np.sqrt(1 + masa_kg / masa_entorno_kg)

    return v_n


### V1: df['V_1']
def modelo_V1(masa_entorno=None, masa_estrella=None, radio=None, v_max=225, describe=None):
    """
    Modelo V1: modelo de materia oscura
    doc: https://es.wikipedia.org/wiki/Perfil_de_Navarro%E2%80%93Frenk%E2%80%93White
    
    Argumentos:
    - masa_entorno: masa total del entorno a la estrella en unidades solares.
    - masa_estrella: masa de la estrella en unidades solares: SOL=1.
    - radio: radio de la órbita en kpc.
    Retorna:
    - Velocidad tangencial en km/s.
    """

    if describe:
        print(modelo_V1.__doc__)
        return
    
    G = 4.3e-6  # Constante gravitacional en unidades kpc^3/(M_sol*Gyr^2)
    masa_total = masa_entorno + masa_estrella
    v_c = np.sqrt(G * masa_total / radio)  # Velocidad circular en km/s
    v_n = v_max * np.sqrt(masa_entorno / masa_total)  # Velocidad de Navarro-Frenk-White en km/s
    v_t = np.sqrt(v_c**2 + v_n**2)  # Velocidad tangencial en km/s

    return v_t


### V2: df['V_2'] Una idea de V. Estrada
def modelo_V2(masa_entorno=None, masa_estrella=None, radio=None, describe=None):
    """
    Modelo V2: Modelo de V. Estrada basado en el modelo einsteniano, 
    que incluye la propagación no instantánea de la fuerza gravitacional.
    El parámetro tiempo_propagacion representa el tiempo de propagación de la fuerza. 
    Este parámetro se utiliza para calcular un factor exponencial 
    que afecta a la velocidad de la estrella en función de la distancia al centro de la galaxia 
    y el tiempo de propagación. El factor exponencial se calcula como 
    math.exp(-radio_al / (v_n * tiempo_propagacion)), 
    donde v_n es la velocidad inicial de la estrella sin el factor correctivo, 
    y radio_al es el radio en años luz.

    Args:
        masa_entorno (float): Masa del entorno o centro de masas de la galaxia en masas solares.
        masa_estrella (float): Masa de la estrella en masas solares.
        radio (float): Distancia desde el centro de masas de la galaxia hasta la estrella en kiloparsecs.
        describe (bool, optional): Si se establece en True, muestra la descripción del modelo. 
                                   Por defecto es None.

    Returns:
        float: Velocidad de la estrella corregida por el factor exponencial.

    ####################################################################
    # La propagación no instantánea de la fuerza gravitacional es una  #
    # idea de Víctor Estrada Díaz que plantea que dadas las distancias #
    # estelares la propagacion de la fuerza ha de jugar un factor de   #
    # vital importancia en las fuerzas y velocidades resultantes.      #
    #                                                                  #
    # Esta idea de Víctor Estrada es un tema abierto de estudio.       #
    #===================================================================

    """
    
    if describe:
        print(modelo_V2.__doc__)
        return
    
    G = 4.3e-6  # Constante gravitacional en unidades kpc^3/(M_sol*Gyr^2)
    masa_total = masa_entorno + masa_estrella
    v_c = np.sqrt(G * masa_total / radio)  # Velocidad circular en km/s
    # ψ''(r) = v^2 * (ψ''(x) + ψ''(y) + ψ''(z)) * ρ(x, y, z)
    v_p=0
    v_t = np.sqrt(v_c**2 + v_p**2)  # Velocidad tangencial en km/s



def modelo_V3(masa_total=None, masa=None, radio=None, describe=None):
    """
    Moselo V3: Modelo original de Einstein con propagación retardada de la fuerza y la deformación del espacio.
    
    Argumentos:
    - masa_total: Masa total de la galaxia en unidades de masa solar.
    - masa: Masa de la estrella en unidades de masa solar.
    - radio: Distancia entre la estrella y el centro galáctico en kiloparsecs (kpc).
    Retorna:
    - Velocidad tangencial considerando la propagación retardada de la fuerza y la deformación del espacio en km/s.
    """

    if describe:
        print(modelo_V3.__doc__)
        return
    
    G = 4.3e-6  # Constante gravitacional en unidades kpc^3/(M_sol*Gyr^2)
    c = 299792.458  # Velocidad de la luz en km/s

    # Cálculo del retardo en la propagación
    retardo = radio * 3.2616  # Conversión de kpc a años luz

    # Cálculo de la fuerza gravitatoria teniendo en cuenta el retardo
    fuerza = G * (masa_total * masa) / (radio**2 + (retardo/c)**2)

    # Cálculo de la velocidad de la estrella
    velocidad = math.sqrt(fuerza / masa)

    return velocidad

def modelo_V4(masa_total=None, masa=None, radio=None, describe=None):
    """
    Modelo V4: Modelo original de Einstein sin velocidad de propagación de la fuerza.
    Este es el primer modelo propuesto por Einstein sin introducir artificios posteriores como la materia oscura.
    
    Fórmulas:
    fuerza = G * (masa_total * masa) / radio**2
    velocidad = math.sqrt(fuerza / masa)
    
    Argumentos:
    - masa_total: Masa total de la galaxia en unidades de masa solar.
    - masa: Masa de la estrella en unidades de masa solar.
    - radio: Distancia entre la estrella y el centro galáctico en kiloparsecs (kpc).
    Retorna:
    - Velocidad tangencial considerando la velocidad de propagación de la fuerza en km/s.
    """

    if describe:
        print(modelo_V3.__doc__)
        return

    G = 4.3e-6  # Constante gravitacional en unidades kpc^3/(M_sol*Gyr^2)

    # Cálculo de la fuerza gravitatoria
    fuerza = G * (masa_total * masa) / radio**2

    # Cálculo de la velocidad de la estrella
    velocidad = math.sqrt(fuerza / masa)

    return velocidad


######### FIN DE MODELOS #########

### funciones auxiliares ###
# MARCA, ponernuna marca en las gráficas
def marca(fig):
    # Definir el texto de la marca
    texto = "© by Victor Estrada Diaz, #1"

    # Obtener las dimensiones de la figura
    ancho_fig = fig.get_figwidth()
    alto_fig = fig.get_figheight()

    # Calcular la posición del texto
    fonts = 8
    x = 0.95
    y = 0.5
    print(f'x:{x},y:{y},ancho:{ancho_fig},alto:{alto_fig}')
    fig.text(x, y, texto, fontsize=fonts, color='lightgray', ha='center', va='center',\
        rotation='vertical',\
        bbox=dict(facecolor='yellow', edgecolor='blue', alpha=0.2), zorder=10)

    
#PLOT función para poner un grafico y asignar un numero de figura
def plot(fig, n):

    f = f'{n}'
    titulo = f'Figura {f}'
    plt.suptitle(titulo)
    marca(fig)
    fig.savefig(f'figura{f}.png')
    plt.clf()
    
#Generar una animación 3d de los datos    
def generate_rotation_animation(df, filename, color):

    if os.path.exists(filename):
        os.remove(filename)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f'Estrellas de la muestra: {df.size}')
    
 # Escalar los valores de masa al rango deseado
    min_masa = df['masa'].min()
    max_masa = df['masa'].max()
    scaled_masa = (df['masa'] - min_masa) / (max_masa - min_masa)  # Escala los valores al rango [0, 1]

    # Asignar los valores escalados a la variable 's' para el tamaño del punto
    s = np.interp(scaled_masa, (0, 1), (0.5, 10))  # Escala el rango [0, 1] a un rango deseado (por ejemplo, de 0.5 a 10)

    scatter = ax.scatter(df['x'], df['y'], df['z'], c=color, s=s)
    ax.scatter(x0, y0, z0, c='r', s=5, zorder=10)


    def progress(i, n):
        print(f'Saving frame {i}/{n}')

    def animate(frame):
        ax.view_init(elev=6, azim=frame * 3)
             
    ani = animation.FuncAnimation(fig, animate, frames=60, interval=50)

    print (f"fichero:{filename}")
    ani.save(filename, writer='pillow')


### INICIO ###

sql_gaia(10000000)

if recalcular or not os.path.exists('estrellas_result.csv'):
    estrellas()
    
    df = leer_estrellas()
    df = completar_estrellas(df,True)

    # Guardar el DataFrame resultante en un archivo csv

    df.to_csv('estrellas_result.csv',header=True)

else:
    df = pd.read_csv('estrellas_result.csv')

    #print(df.iloc[0:6, -5:])
    modelo_V1(describe=True)
    modelo_V2(describe=True)
    modelo_V3(describe=True)
    modelo_V4(describe=True)
    
    df = completar_estrellas(df, False)

### PRESENTAR RESULTADOS ###

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#definir que graficas trataremos
todas = True #False #Si true las generará todas
nfig = 4      #poner la grafica a trabajar
if todas:
    nfig = 0

#borra los archivos de imagen, solo si graficamos todas, si no la que trabajemos
archivos = []
for i in range(0,20):
    if todas or nfig == i:
        archivos = archivos + [f"figura{i}.png"]
for archivo in archivos:
    try:
        print ("borrando ... ",archivo)
        os.remove(archivo)
    except FileNotFoundError:
        continue
        

### Color de las estrellas ###
print("MUestras en la galaxia: XYZ\n")
print(" Color indica tamaño de la estrella\n")

# Definir el tamaño de la cuadrícula
grid_size = 100

# Calcular las coordenadas de la cuadrícula
x_min, x_max = df['x'].min(), df['x'].max()
y_min, y_max = df['y'].min(), df['y'].max()
z_min, z_max = df['z'].min(), df['z'].max()
x_grid = np.linspace(x_min, x_max, grid_size)
y_grid = np.linspace(y_min, y_max, grid_size)
z_grid = np.linspace(z_min, z_max, grid_size)

# Calcular la densidad de puntos
density, _ = np.histogramdd([df['x'], df['y'], df['z']], bins=[x_grid, y_grid, z_grid])
print (f"Estrellas:{df['x'].size} --- histograma:{density.size}")

# Normalizar la densidad en el rango 0-1
density_norm = (density - density.min()) / (density.max() - density.min())

# Colores de la nuve de puntos
# Transformar las masas a escala logarítmica
masas_log = np.log10(df['masa'])
vmin = np.log10(df['masa'].min())
vmax = np.log10(df['masa'].max())
# Calcular los colores en la escala logarítmica
colors = cm.viridis((masas_log - vmin) / (vmax - vmin))

#### Gráfica ####
#### Estrellas muestreadas de la galaxia.
# Crear el gráfico 3D
# Generar y guardar la animación combinada
if not os.path.exists('estrellas.gif'): generate_rotation_animation(df, 'estrellas.gif',colors)


if todas or nfig==0:
    ### Gráfica ###
    ### Velocidad tangencial medida desde las observaciones
    fig = plt.figure()
    ax = fig.add_subplot(111)
    rango_y=(0,1200) #Km/seg 
    ax.set_ylim(rango_y)
    ax.scatter(df['radio'], df['v_6'],c=colors,s=0.1)
    ax.set_ylabel('V6 (Km/seg)')
    ax.set_xlabel('radio')
    ax.set_title(f'Velocidad tangencial real GSM desde las observaciones de Gaia')         
    ax.scatter(d0, v0, c='r', marker='o', zorder=10)
    plot(fig,nfig)
    if todas: nfig += 1

if todas or nfig==1:
    ### Gráfica ### Histograma
    ### Distribución de las masas vistas
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # plotear el histograma 
    sns.kdeplot(df['masa'], ax=ax, color='blue', linewidth=2, alpha=0.5)
    ax.set_ylabel('distribucion')
    ax.set_xlabel('masa')
    ax.set_title(f'Masa estelar y distribución')         
    plot(fig,nfig)
    if todas: nfig += 1

if todas or nfig==2:
    ### Gráfica ###
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(df['radio'], df['masa'],c=colors,s=0.1)
    ax.set_ylabel('masa')
    ax.set_xlabel('radio')
    ax.set_title(f'Masa estelar y localización')         
    ax.scatter(d0, 1, c='r', marker='o')
    plot(fig,nfig)
    if todas: nfig += 1
    
if todas or nfig==3:
    #### Grafica ### 
    #### Velocidad según la masa
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Ajustar los límites del eje y
    rango_y=(0,2400) #Km/seg 
    ax.set_ylim(rango_y)
    ax.scatter(df['masa'], df['v_6'],c='purple',s=0.1, zorder=5)
    ax.set_ylabel('V6')
    ax.set_xlabel('masa')
    ax.set_title(f'Velocidad tangencial según masa')         
    ax.scatter([d0], [v0], c='r', marker='o', s=10)
    plot(fig,nfig)
    if todas: nfig += 1

if todas or nfig==4:
    ### Gráfica ###
    ### Velocidad tangencial medida desde las observaciones
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #rango_y=(0,1200) #Km/seg 
    #ax.set_ylim(rango_y)
    #ax.scatter(df['radio'], df['v_6'],c=colors,s=0.1,zorder=1)
    #ax.scatter(df['radio'], df['v_7'],c='black',s=0.1, zorder=9)
    ax.scatter(df['radio'], df['v_g'],c='black',s=0.1, zorder=9)
    ax.set_ylabel('V6, Vg (Km/seg)')
    ax.set_xlabel('radio')
    ax.set_title(f'Velocidad tangencial contra Modelo 4 Mejor ajuste')         
    #ax.scatter(d0, v0, c='r', marker='o', zorder=10)
    plot(fig,nfig)
    if todas: nfig += 1

if todas or nfig==5:
    ### Gráfica ###
    ### Velocidad tangencial por modelo Ve
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Ajustar los límites del eje y
    rango_y=(0,1200) #Km/seg 
    ax.set_ylim(rango_y)
    ax.scatter(df['radio'], df['v_6'],c=colors,s=0.1)
    ax.scatter(df['radio'], df['v_4'],c='y',s=0.1)
    ax.set_ylabel('V6, V4 (Km/seg)')
    ax.set_xlabel('radio')
    ax.set_title(f'Velocidad tangencial contra modelo V4 de Einstein sin materia oscura')
    ax.scatter(d0, v0, c='r', marker='o')
    plot(fig,nfig)
    if todas: nfig += 1

if todas or nfig==6:
    ### Gráfica ###
    ### Velocidad tangencial por modelo V1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Ajustar los límites del eje y
    rango_y=(0,1200) #Km/seg 
    ax.set_ylim(rango_y)
    ax.scatter(df['radio'], df['v_6'],c=colors,s=0.1)
    ax.scatter(df['radio'], df['v_1'],c='y',s=0.1)
    ax.set_ylabel('V_6, V1 (Km/seg)')
    ax.set_xlabel('radio')
    ax.set_title(f'Velocidad tangencial V6 contra el modelo V1')         
    ax.scatter(d0, v0, c='r', marker='o')
    plot(fig,nfig)
    if todas: nfig += 1

if todas or nfig==7:
    ### Gráfica ###
    ### Velocidad tangencial por modelo V1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Ajustar los límites del eje y
    rango_y=(0,1200) #Km/seg 
    ax.set_ylim(rango_y)
    ax.scatter(df['radio'], df['v_6'],c=colors,s=0.1)
    ax.scatter(df['radio'], df['v_n'],c='y',s=0.1)
    ax.set_ylabel('VN')
    ax.set_xlabel('radio')
    ax.set_title(f'Velocidad tangencial V6 contra modelo Vn de Newton')         
    ax.scatter(d0, v0, c='r', marker='o')
    plot(fig,nfig)
    if todas: nfig += 1

if todas or nfig==8:
    ### Gráfica ###
    ### Velocidad tangencial por modelo V1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Ajustar los límites del eje y
    rango_y=(0,1200) #Km/seg 
    ax.set_ylim(rango_y)
    ax.scatter(df['radio'], df['v_6'],c=colors,s=0.1)
    ax.scatter(df['radio'], df['v_1'],c='b',s=0.1)
    ax.set_ylabel('V6, V1 (Km/seg)')
    ax.set_xlabel('radio')
    ax.set_title(f'Velocidad tangencial real v6 contra modelo V1')         
    ax.scatter(d0, v0, c='r', marker='o')
    plot(fig,nfig)
    if todas: nfig += 1

if todas or nfig==9:
    ### Gráfica ###
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Ajustar los límites del eje y
    rango_y=(0,1200) #Km/seg 
    ax.set_ylim(rango_y)
    ax.scatter(df['radio'], df['v_2'],c='purple',s=0.1, zorder=10)
    ax.scatter(df['radio'], df['vt'],c=colors,s=0.1)
    ax.set_ylabel('V2 vs VT')
    ax.set_xlabel('radio')
    ax.set_title(f'Velocidad real contra modelo V2 propaga fuerza')         
    ax.scatter([d0], [v0], c='r', marker='o', s=10)
    plot(fig,nfig)
    if todas: nfig += 1

if todas or nfig==10:
    ### Gráfica  
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Ajustar los límites del eje y
    rango_y=(0,1200) #Km/seg 
    ax.set_ylim(rango_y)
    ax.scatter(df['radio'], df['v_n'],c='yellow',s=0.1, zorder=9)
    ax.scatter(df['radio'], df['v_1'],c='orange',s=0.1, zorder=9)
    ax.scatter(df['radio'], df['v_2'],c='black',s=0.1, zorder=9)
    ax.scatter(df['radio'], df['v_e'],c='blue',s=0.1, zorder=9)
    ax.scatter(df['radio'], df['v_3'],c='purple',s=0.1, zorder=9)
    ax.scatter(df['radio'], df['v_4'],c='cyan',s=0.1, zorder=9)
    ax.scatter(df['radio'], df['v_5'],c='red',s=0.1, zorder=9)
    ax.set_ylabel('v1:o, v2:n, Ve:b, V3:p, V4:c')
    ax.set_xlabel('radio')
    ax.set_title(f'Velocidad calculad por los modelos por estrella')         
    ax.scatter([d0], [v0], c='r', marker='o', s=10, zorder=10)
    plot(fig,nfig)
    if todas: nfig += 1

if todas or nfig==11:
    ### Gráfica ###
    ### Velocidad tangencial por mi modelo V2 contra real VT
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Ajustar los límites del eje y
    rango_y=(0,1200) #Km/seg 
    ax.set_ylim(rango_y)
    ax.scatter(df['radio'], df['v_n'],c='yellow',s=0.1, zorder=10)
    ax.scatter(df['radio'], df['v_1'],c='orange',s=0.1, zorder=10)
    ax.scatter(df['radio'], df['v_2'],c='black',s=0.1, zorder=10)
    ax.scatter(df['radio'], df['v_e'],c='blue',s=0.1, zorder=10)
    ax.scatter(df['radio'], df['v_3'],c='purple',s=0.1, zorder=10)
    ax.scatter(df['radio'], df['vt'],c=colors,s=0.1)
    ax.scatter(df['radio'], df['v_4'],c='cyan',s=0.1, zorder=10)
    ax.set_ylabel('Vn:y, v1:o, v2:n, Ve:b, V3:p, V4:c Km/seg')
    ax.set_xlabel('radio')
    ax.set_title(f'Velocidad tangencial real comparada con modelos')         
    ax.scatter([d0], [v0], c='r', marker='o', s=10)
    plot(fig,nfig)
    fig = plt.figure()
    if todas: nfig += 1
    
if todas or nfig==12:
    ###Grafica###
    ax = fig.add_subplot(111)
    # Ajustar los límites del eje y
    rango_y=(0,1200) #Km/seg 
    ax.set_ylim(rango_y)
    ax.scatter(df['radio'], df['v_6'],s=0.1, zorder=8,c=colors)
    ax.scatter(df['radio'], df['v_1'],c='orange',s=0.1, zorder=9)
    ax.scatter(df['radio'], df['v_2'],c='green',s=0.1, zorder=10)
    ax.scatter(df['radio'], df['v_3'],c='blue',s=0.1, zorder=7)
    ax.scatter(df['radio'], df['v_4'],c='cyan',s=0.1, zorder=9)
    #ax.scatter(df['radio'], df['v_5'],c='cyan',s=0.1, zorder=9)
    ax.scatter(df['radio'], df['v_n'],c='yellow',s=0.1, zorder=9)
    ax.set_ylabel('vGSM:color, v1:ora, v2:gre, v3:yel, Km/s')
    ax.set_xlabel('radio')
    ax.set_title(f'Velocidad GSM')         
    ax.scatter([d0], [v0], c='r', marker='o', s=10, zorder=10)
    plot(fig,nfig)
    if todas: nfig += 1
