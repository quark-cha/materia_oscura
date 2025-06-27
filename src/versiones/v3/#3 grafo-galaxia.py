pip import sys
from pathlib import Path

# AÑADIR LA RAÍZ DEL PROYECTO A sys.path
ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))


import pandas as pd
import networkx as nx
import pickle
import numpy as np
import matplotlib.pyplot as plt
import math
import shutil
import ast

# 1. Importar el sistema de verificación
from dependencias import verificar_dependencias

# 2. Verificar dependencias al inicio del programa
problemas = verificar_dependencias()
import sys
if problemas:
        
        sys.exit(1)
        

# AÑADIR LA RAÍZ DEL PROYECTO A sys.path
from pathlib import Path 
ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))


import pandas as pd
import networkx as nx
import pickle
import numpy as np
import matplotlib.pyplot as plt
import math
import shutil
import ast


#otros modulos del proyecto: en _VED_.PY
global grafo
grafo = nx.Graph()

from src._ved_ import ver, plot

class Objeto:
    def __init__(self):
        self.tipo = ""

_ved_ = Objeto()
_ved_.tipo = "galaxia"


archivo_grafo = 'datos/grafo_estrellas.gpickle'
archivo_copia = 'datos/grafo_estrellas_backup.gpickle'
archivo_grafo_csv = 'datos/grafo.csv'


def vector_v(ra, dec, parallax, pm_ra, pm_dec, radial_velocity, test=False):
    """_summary_

    Args:
        ra (_type_): ascension recta
        dec (_type_): declinacion
        parallax (_type_): paralelage
        pm_ra (_type_): movimiento aparente en (ra)
        pm_dec (_type_): movimiento aparente en (dec)
        radial_velocity (_type_): velocidad radial
        test (bool, optional): _description_. Defaults to False.

    Returns:
        'vx': [vx_real], velocidad real GSM en la direccion X
        'vy': [vy_real], velocidad real GSM en la direccion Y
        'vz': [vz_real], velocidad real GSM en la direccion Z
        'mod_v': [mod_v], modulo de la velocidad real XYZ
        'mod_t': [mod_t], modulo de la velocidad real enla tangente
        'divergencia': [divergencia] diferencia de modulos
    """
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
        'vx': vx_real,
        'vy': vy_real,
        'vz': vz_real,
        'mod_v': mod_v,
        'mod_t': mod_t,
        'divergencia': divergencia
    }
    
    return data

def evolucionar_grafo(G,recalcular_velocidades= False):
    #comprueba si el grafo necesita datos que faltan o precisa evolucionar
    
    print ("¿evolucionar_grafo?")
    porcentaje_completado = 0.0
    contador=0
    num_estrellas=G.number_of_nodes()
    n = list(G.nodes)[0]  # Accede al primer nodo en el grafo
    if 'mod_v' not in G.nodes[n] or recalcular_velocidades:
        for estrella in G.nodes:
            if recalcular_velocidades:
                # Eliminar las columnas existentes antes de recalcular
                cols_to_remove = ['vx', 'vy', 'vz', 'mod_v', 'mod_t', 'divergencia']
                for col in cols_to_remove:
                    if col in G.nodes[estrella]:
                        G.nodes[estrella].pop(col)
            masa = G.nodes[estrella]['masa']
            x= G.nodes[estrella]['x']
            y= G.nodes[estrella]['y']
            z= G.nodes[estrella]['z']
            radio = G.nodes[estrella]['radio']
            radial_velocity = G.nodes[estrella]['vv']
            ra = G.nodes[estrella]['ra']
            dec = G.nodes[estrella]['dec']
            parallax = G.nodes[estrella]['parallax']
            pmra = G.nodes[estrella]['pmra']
            pmdec = G.nodes[estrella]['pmdec']
            v = vector_v(ra, dec, parallax, pmra, pmdec, radial_velocity)
            G.nodes[estrella]['vx']=v['vx']
            G.nodes[estrella]['vy']=v['vy']
            G.nodes[estrella]['vz']=v['vz']
            G.nodes[estrella]['mod_v']=v['mod_v']
            G.nodes[estrella]['mod_t']=v['mod_t']
            G.nodes[estrella]['divergencia']=v['divergencia']
            
            contador += 1  # Incrementar el contador

            # Calcular el porcentaje completado
            test = porcentaje_completado
            porcentaje_completado = ((contador / num_estrellas) * 100) // 1

            # Imprimir el porcentaje en la misma línea
            if test != porcentaje_completado: 
                print(f"Progreso de evolucion: {porcentaje_completado:.2f}%")

    return G

# 3. Guardar el nuevo grafo con los atributos adicionales
import pandas as pd

def guardar_grafo_a_csv(G, archivo_csv):
    # Crear una lista de nodos con sus atributos
    nodo_atributos = [(nodo, G.nodes[nodo]) for nodo in G.nodes]
    
    # Convertir la lista de nodos y atributos en un DataFrame
    df = pd.DataFrame(nodo_atributos)
    
    # Guardar el DataFrame en un archivo CSV
    df.to_csv(archivo_csv, index=False, sep=';')
    
    print(f"Guardado grafo en {archivo_csv} con estrellas {len(G.nodes)}")


import ast

def cargar_grafo_desde_csv(grafo_csv):
    print("Creando un nuevo grafo desde ", grafo_csv)
    df = pd.read_csv(grafo_csv)
    # Aplicar filtros para eliminar datos con masa nula o mayor a 10000 soles
    df_filtrado = df.dropna(subset=['masa'])  # Eliminar filas con masa nula
    df_filtrado = df_filtrado[df_filtrado['masa'] <= 10000]  # Filtrar valores mayores a 10000
    df=df_filtrado

    G = nx.Graph()

    for _, row in df.iterrows():
        etiqueta = row[df.columns[0]]
        datos_nodo = {col: ast.literal_eval(row[col]) for col in df.columns[1:]}
        G.add_node(etiqueta, **datos_nodo)

    print(f"Creado grafo en {grafo_csv} con {G.number_of_nodes()} estrellas")
    return G


def crear_grafo():
    # Crea un grafo totalmente nuevo desde un archivo csv con datos iniciales
    # Leer el archivo CSV con los datos de las estrellas generado desde GAIA
    df = pd.read_csv('datos/estrellas.csv')
    df_filtrado = df.dropna(subset=['masa'])  # Eliminar filas con masa nula
    df_filtrado = df_filtrado[df_filtrado['masa'] <= 2000]  
    # Filtrar valores mayores a 2000 hay valores de 2.88e+03 
    data=df_filtrado
    del df_filtrado
    del df
    
    print(f'Leidas {data.size} estrellas de datos/estrellas.csv')
    print("Ahora vamos a hacer unos calculos preliminares")

    # Crear un grafo vacío
    grafo = nx.Graph()

    # Agregar nodos al grafo
    for index, row in data.iterrows():
        import numpy as np
        masa = row['masa']  # Obtener la masa de la estrella del dataframe
        
        # Agregar el nodo al grafo con la masa conocida

        nodo_id = index  # Puedes usar el índice como identificador del nodo
        x= row['x']
        y= row['y']
        z= row['z']
        
        radio = row['radio']
        radial_velocity = row['vv']
        ra = row['ra']
        dec = row['dec']
        parallax = row['parallax']
        pmra = row['pmra']
        pmdec = row['pmdec']
        v = vector_v(ra, dec, parallax, pmra, pmdec, radial_velocity)

        # Agregar atributos adicionales al nodo
        atributos = {
            'x': x,                 #distancia X al centro galáctico
            'y': y,                 #distancia Y al centro galáctico
            'z': z,                 #distancia Z al centro galáctico
            'masa': masa,           #masa en soles de la estrella
            'radio': radio,         #distancia al centro galactico o radio
            'vv': radial_velocity,  #velocidad radial estrimada por Gaia
            'ra': ra,               #coordenada en ascensión recta
            'dec': dec,             #coordenada en declinacion
            'parallax': parallax,   #coordenada en paralelaje
            'pmra': pmra,           #movimiento aparente en ascension recta
            'pmdec'
            
            : pmdec,         #movimiento aparente en declinación
            'vx': v['vx'],              #velocidad real (GSM) en coordenada X
            'vy': v['vy'],              #velocidad real (GSM) en coordenada Y
            'vz': v['vz'],              #velocidad real (GSM) en coordenada Z
            'entorno': 0,           #numero de estrellas en el entorno de influencia de la estrella
            'dif_vx':0.0,           #diferencia velocidad real - velocidad teorica t1 coordenada X
            'dif_vy':0.0,           #diferencia velocidad real - velocidad teorica t1 coordenada Y
            'dif_vz':0.0,           #diferencia velocidad real - velocidad teorica t1 coordenada Z
            'vx_t1': 0.0,           #velocidad teorica t1 en coordenada X
            'vy_t1': 0.0,           #velocidad teorica t1 en coordenada Y
            'vz_t1': 0.0,           #velocidad teorica t1 en coordenada Z
            'mod_v': v['mod_v'],           #modulo de la velocidad real (GSM)
            'mod_t': v['mod_t'],           #modulo de la velocidad real (GSM) tangencial
            'divergencia':v['mod_v']-v['mod_t']   #divergencia velocidad teorica t1 - velocidad real
            }
     
        grafo.add_node(nodo_id, **atributos)

    # Agregar aristas al grafo (conexiones entre las estrellas)
    # ... Código para agregar aristas según tu criterio ...

    print ("Grafo Creado")

    return grafo

def actualizar_grafo_en_disco(grafo, archivo):
    with open(archivo, 'wb') as f:
        pickle.dump(grafo, f)

    
# Constante de gravitación universal (m^3/kg/s^2)
G = 6.67430e-11
# Masa de la estrella grande en kilogramos
masa_grande = 1000 * 1.9885e30  # 1000 masas solares en kg
umbral_velocidad = 1.0  # m/seg Umbral de velocidad
kpc_a_metros = 3.086e19  # Kiloparsec a metros

# Calcular la distancia crítica en kiloparsecs:
# r = sqrt(G * M^2 * t^2 / (v * M)) = sqrt(G * M / v) en metros
#distancia_critica_metros = math.sqrt(G * masa_grande / umbral_velocidad) #
distancia_critica_kpc = 0.2 # 652,313 años luz

global exc
exc=0
def calcular_fuerza(masa1, masa2, distancia, einstein=False):
    """
    Calcula la fuerza gravitatoria entre dos estrellas utilizando la ley de gravitación de Newton o
    una aproximación simplificada basada en la relatividad general.

    Parámetros:
        masa1 (float): Masas solares.
        masa2 (float): Masas solares.
        distancia (float): En kpc.
        einstein (bool, opcional): Si es True, utiliza la aproximación basada en la relatividad general.

    Retorna:
        float: Fuerza gravitatoria entre las dos estrellas en newtons.
    """
    # Constante de gravitación universal (m^3/kg/s^2)
    constante_gravitacion = 6.67430e-11
    # Velocidad de la luz (m/s)
    velocidad_luz = 299792458.0
    
    #convertir unidades
    masa1_kg = masa1 * 1.9885e30  # Convertir masas solares a kg
    masa2_kg = masa2 * 1.9885e30
    distancia_metros = distancia * 3.086e19  # kpc a metros
    
    global exc
    
    if einstein:   
        # Fórmula basada en la relatividad general para la "fuerza equivalente"
        fuerza = (masa1_kg * masa2_kg * constante_gravitacion) / (distancia_metros**2 * velocidad_luz**2)
    else:
        # Fórmula clásica de la ley de gravitación de Newton
        fuerza = 0.0
        try:
            fuerza = constante_gravitacion * (1 / distancia_metros**2) * masa2_kg * masa1_kg
        except:
            exc += 1
            print(f"excepcion [{exc}] m1:{masa1:.2e} m2:{masa2:.2e}")
        # fuerza en Newtons

    return fuerza


def velocidad(grafo):
    estrellas = grafo.nodes
    num_estrellas = len(estrellas)

    # Constantes
    masa_del_sol_kg = 1.9885e30  # Masa del Sol en kilogramos
    kpc_a_metros = 3.086e19  # Kiloparsec a metros

    # Iterar a través de las estrellas
    contador = 0
    porcentaje_completado = 0
        
    for i, estrella in enumerate(estrellas):
        for j, estrella_grande in enumerate(estrellas):
            if j >= i:  # Solo calcular para las estrellas por encima de la diagonal
                continue
                    # Calcula la distancia entre las estrellas estrella1 y estrella2
            distancia_kpc = np.sqrt \
                        (   (estrellas[estrella]['x'] - estrellas[estrella_grande]['x'])**2 + 
                            (estrellas[estrella]['y'] - estrellas[estrella_grande]['y'])**2 + 
                            (estrellas[estrella]['z'] - estrellas[estrella_grande]['z'])**2)

            # Comprobar si la distancia cumple con el criterio
            if distancia_kpc > distancia_critica_kpc:
                continue
            
            estrellas[estrella]['entorno'] += 1
            estrellas[estrella_grande]['entorno'] += 1
                        
            Masa_sol = masa_del_sol_kg
            M1 = estrellas[estrella]['masa']
            M2 = estrellas[estrella_grande]['masa']
            c_kpc_a_metros = kpc_a_metros
            fuerza = calcular_fuerza(M1, M2, distancia_kpc)  # Fuerza en newtons
            m1_x = estrellas[estrella]['x']*kpc_a_metros
            m1_y = estrellas[estrella]['y']*kpc_a_metros
            m1_z = estrellas[estrella]['z']*kpc_a_metros
            m2_x = estrellas[estrella_grande]['x']*kpc_a_metros
            m2_y = estrellas[estrella_grande]['y']*kpc_a_metros
            m2_z = estrellas[estrella_grande]['z']*kpc_a_metros
            
            # Calcular componentes de la fuerza en x, y, z
            fuerza_x = fuerza * abs(m1_x - m2_x) / (c_kpc_a_metros * distancia_kpc) 
            fuerza_y = fuerza * abs(m1_y - m2_y) / (c_kpc_a_metros * distancia_kpc) 
            fuerza_z = fuerza * abs(m1_z - m2_z) / (c_kpc_a_metros * distancia_kpc)

            # La velocidad que aporta es la fuerza dividida por su masa
            # y pasamos estas velocidades de m/seg a km/seg
            # Actualizar velocidades de las estrellas
            estrellas[estrella]['vx_t1'] += (fuerza_x / (M1*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s
            estrellas[estrella]['vy_t1'] += (fuerza_y / (M1*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s
            estrellas[estrella]['vz_t1'] += (fuerza_z / (M1*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s

            # Aplicar la acción inversa en la estrella grande
            estrellas[estrella_grande]['vx_t1'] -= (fuerza_x / (M2*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s
            estrellas[estrella_grande]['vy_t1'] -= (fuerza_y / (M2*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s
            estrellas[estrella_grande]['vz_t1'] -= (fuerza_z / (M2*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s            
        
        
        
        contador += 1  # Incrementar el contador

        # Calcular el porcentaje completado
        test = porcentaje_completado
        porcentaje_completado = ((contador / num_estrellas) * 100) // 1

        # Imprimir el porcentaje en la misma línea
        if test != porcentaje_completado: 
            print(f"Progreso: {porcentaje_completado:.2f}%, ultima entorno:{estrellas[estrella]['entorno']}")

        # Esperar un momento para simular el progreso
        # time.sleep(0.1)  # Puedes ajustar este valor para controlar la velocidad de actualización

    # Calcular aceleraciones y actualizar velocidades para cada estrella
    for i, estrella in enumerate(estrellas):
        # Diferencia entre la velocidad teórica y la velocidad efectiva (GSM)
        estrellas[estrella]['dif_vx'] = estrellas[estrella]['vx'] - estrellas[estrella]['vx_t1']
        estrellas[estrella]['dif_vy'] = estrellas[estrella]['vy'] - estrellas[estrella]['vy_t1']
        estrellas[estrella]['dif_vz'] = estrellas[estrella]['vz'] - estrellas[estrella]['vz_t1']
        
    guardar_grafo_a_csv(grafo, grafo_csv)

    # Imprimir mensaje de finalización
    print("\nCálculo completado.")




archivo_grafo = 'datos/grafo_estrellas.gpickle'
archivo_copia = 'datos/grafo_estrellas_backup.gpickle'
grafo_csv = 'datos/grafo.csv'

# Comprobar si el archivo del grafo existe
try:
    with open(archivo_grafo, 'rb') as f:
        grafo = pickle.load(f)
        print(f"1 - procesado binario {archivo_grafo}, {grafo.number_of_nodes()} estrellas")
        #comprobar si faltan calculos y completarlos
        grafo=evolucionar_grafo(grafo, True)     
        guardar_grafo_a_csv(grafo, grafo_csv)
                               
    # Crear una copia de seguridad del archivo del grafo
    shutil.copyfile(archivo_grafo, archivo_copia)
     
except FileNotFoundError:
    
    try:
        grafo=cargar_grafo_desde_csv(grafo_csv)
        print(f"2 - procesado grafo_csv {grafo_csv}, {grafo.number_of_nodes()} estrellas")        
        
    except FileNotFoundError:
        # si nada es posible entoces crea un grafo nuevo a partir de los datos de Gaia
        print("El archivo del grafo no existe. Creando el grafo...")
        grafo = crear_grafo()
        velocidad(grafo)
        print(f"3 - creado nuevo, {grafo.number_of_nodes()} estrellas")                

actualizar_grafo_en_disco(grafo, archivo_grafo)

ver(grafo)
print("Trabajo terminado")


#otros modulos del proyecto: en _VED_.PY
global grafo
grafo = nx.Graph()

from src._ved_ import ver, plot

class Objeto:
    def __init__(self):
        self.tipo = ""

_ved_ = Objeto()
_ved_.tipo = "galaxia"


archivo_grafo = 'datos/grafo_estrellas.gpickle'
archivo_copia = 'datos/grafo_estrellas_backup.gpickle'
archivo_grafo_csv = 'datos/grafo.csv'


def vector_v(ra, dec, parallax, pm_ra, pm_dec, radial_velocity, test=False):
    """_summary_

    Args:
        ra (_type_): ascension recta
        dec (_type_): declinacion
        parallax (_type_): paralelage
        pm_ra (_type_): movimiento aparente en (ra)
        pm_dec (_type_): movimiento aparente en (dec)
        radial_velocity (_type_): velocidad radial
        test (bool, optional): _description_. Defaults to False.

    Returns:
        'vx': [vx_real], velocidad real GSM en la direccion X
        'vy': [vy_real], velocidad real GSM en la direccion Y
        'vz': [vz_real], velocidad real GSM en la direccion Z
        'mod_v': [mod_v], modulo de la velocidad real XYZ
        'mod_t': [mod_t], modulo de la velocidad real enla tangente
        'divergencia': [divergencia] diferencia de modulos
    """
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
        'vx': vx_real,
        'vy': vy_real,
        'vz': vz_real,
        'mod_v': mod_v,
        'mod_t': mod_t,
        'divergencia': divergencia
    }
    
    return data

def evolucionar_grafo(G,recalcular_velocidades= False):
    #comprueba si el grafo necesita datos que faltan o precisa evolucionar
    
    print ("¿evolucionar_grafo?")
    porcentaje_completado = 0.0
    contador=0
    num_estrellas=G.number_of_nodes()
    n = list(G.nodes)[0]  # Accede al primer nodo en el grafo
    if 'mod_v' not in G.nodes[n] or recalcular_velocidades:
        for estrella in G.nodes:
            if recalcular_velocidades:
                # Eliminar las columnas existentes antes de recalcular
                cols_to_remove = ['vx', 'vy', 'vz', 'mod_v', 'mod_t', 'divergencia']
                for col in cols_to_remove:
                    if col in G.nodes[estrella]:
                        G.nodes[estrella].pop(col)
            masa = G.nodes[estrella]['masa']
            x= G.nodes[estrella]['x']
            y= G.nodes[estrella]['y']
            z= G.nodes[estrella]['z']
            radio = G.nodes[estrella]['radio']
            radial_velocity = G.nodes[estrella]['vv']
            ra = G.nodes[estrella]['ra']
            dec = G.nodes[estrella]['dec']
            parallax = G.nodes[estrella]['parallax']
            pmra = G.nodes[estrella]['pmra']
            pmdec = G.nodes[estrella]['pmdec']
            v = vector_v(ra, dec, parallax, pmra, pmdec, radial_velocity)
            G.nodes[estrella]['vx']=v['vx']
            G.nodes[estrella]['vy']=v['vy']
            G.nodes[estrella]['vz']=v['vz']
            G.nodes[estrella]['mod_v']=v['mod_v']
            G.nodes[estrella]['mod_t']=v['mod_t']
            G.nodes[estrella]['divergencia']=v['divergencia']
            
            contador += 1  # Incrementar el contador

            # Calcular el porcentaje completado
            test = porcentaje_completado
            porcentaje_completado = ((contador / num_estrellas) * 100) // 1

            # Imprimir el porcentaje en la misma línea
            if test != porcentaje_completado: 
                print(f"Progreso de evolucion: {porcentaje_completado:.2f}%")

    return G

# 3. Guardar el nuevo grafo con los atributos adicionales
import pandas as pd

def guardar_grafo_a_csv(G, archivo_csv):
    # Crear una lista de nodos con sus atributos
    nodo_atributos = [(nodo, G.nodes[nodo]) for nodo in G.nodes]
    
    # Convertir la lista de nodos y atributos en un DataFrame
    df = pd.DataFrame(nodo_atributos)
    
    # Guardar el DataFrame en un archivo CSV
    df.to_csv(archivo_csv, index=False, sep=';')
    
    print(f"Guardado grafo en {archivo_csv} con estrellas {len(G.nodes)}")


import ast

def cargar_grafo_desde_csv(grafo_csv):
    print("Creando un nuevo grafo desde ", grafo_csv)
    df = pd.read_csv(grafo_csv)
    # Aplicar filtros para eliminar datos con masa nula o mayor a 10000 soles
    df_filtrado = df.dropna(subset=['masa'])  # Eliminar filas con masa nula
    df_filtrado = df_filtrado[df_filtrado['masa'] <= 10000]  # Filtrar valores mayores a 10000
    df=df_filtrado

    G = nx.Graph()

    for _, row in df.iterrows():
        etiqueta = row[df.columns[0]]
        datos_nodo = {col: ast.literal_eval(row[col]) for col in df.columns[1:]}
        G.add_node(etiqueta, **datos_nodo)

    print(f"Creado grafo en {grafo_csv} con {G.number_of_nodes()} estrellas")
    return G


def crear_grafo():
    # Crea un grafo totalmente nuevo desde un archivo csv con datos iniciales
    # Leer el archivo CSV con los datos de las estrellas generado desde GAIA
    df = pd.read_csv('datos/estrellas.csv')
    df_filtrado = df.dropna(subset=['masa'])  # Eliminar filas con masa nula
    df_filtrado = df_filtrado[df_filtrado['masa'] <= 2000]  
    # Filtrar valores mayores a 2000 hay valores de 2.88e+03 
    data=df_filtrado
    del df_filtrado
    del df
    
    print(f'Leidas {data.size} estrellas de datos/estrellas.csv')
    print("Ahora vamos a hacer unos calculos preliminares")

    # Crear un grafo vacío
    grafo = nx.Graph()

    # Agregar nodos al grafo
    for index, row in data.iterrows():
        import numpy as np
        masa = row['masa']  # Obtener la masa de la estrella del dataframe
        
        # Agregar el nodo al grafo con la masa conocida

        nodo_id = index  # Puedes usar el índice como identificador del nodo
        x= row['x']
        y= row['y']
        z= row['z']
        
        radio = row['radio']
        radial_velocity = row['vv']
        ra = row['ra']
        dec = row['dec']
        parallax = row['parallax']
        pmra = row['pmra']
        pmdec = row['pmdec']
        v = vector_v(ra, dec, parallax, pmra, pmdec, radial_velocity)

        # Agregar atributos adicionales al nodo
        atributos = {
            'x': x,                 #distancia X al centro galáctico
            'y': y,                 #distancia Y al centro galáctico
            'z': z,                 #distancia Z al centro galáctico
            'masa': masa,           #masa en soles de la estrella
            'radio': radio,         #distancia al centro galactico o radio
            'vv': radial_velocity,  #velocidad radial estrimada por Gaia
            'ra': ra,               #coordenada en ascensión recta
            'dec': dec,             #coordenada en declinacion
            'parallax': parallax,   #coordenada en paralelaje
            'pmra': pmra,           #movimiento aparente en ascension recta
            'pmdec'
            
            : pmdec,         #movimiento aparente en declinación
            'vx': v['vx'],              #velocidad real (GSM) en coordenada X
            'vy': v['vy'],              #velocidad real (GSM) en coordenada Y
            'vz': v['vz'],              #velocidad real (GSM) en coordenada Z
            'entorno': 0,           #numero de estrellas en el entorno de influencia de la estrella
            'dif_vx':0.0,           #diferencia velocidad real - velocidad teorica t1 coordenada X
            'dif_vy':0.0,           #diferencia velocidad real - velocidad teorica t1 coordenada Y
            'dif_vz':0.0,           #diferencia velocidad real - velocidad teorica t1 coordenada Z
            'vx_t1': 0.0,           #velocidad teorica t1 en coordenada X
            'vy_t1': 0.0,           #velocidad teorica t1 en coordenada Y
            'vz_t1': 0.0,           #velocidad teorica t1 en coordenada Z
            'mod_v': v['mod_v'],           #modulo de la velocidad real (GSM)
            'mod_t': v['mod_t'],           #modulo de la velocidad real (GSM) tangencial
            'divergencia':v['mod_v']-v['mod_t']   #divergencia velocidad teorica t1 - velocidad real
            }
     
        grafo.add_node(nodo_id, **atributos)

    # Agregar aristas al grafo (conexiones entre las estrellas)
    # ... Código para agregar aristas según tu criterio ...

    print ("Grafo Creado")

    return grafo

def actualizar_grafo_en_disco(grafo, archivo):
    with open(archivo, 'wb') as f:
        pickle.dump(grafo, f)

    
# Constante de gravitación universal (m^3/kg/s^2)
G = 6.67430e-11
# Masa de la estrella grande en kilogramos
masa_grande = 1000 * 1.9885e30  # 1000 masas solares en kg
umbral_velocidad = 1.0  # m/seg Umbral de velocidad
kpc_a_metros = 3.086e19  # Kiloparsec a metros

# Calcular la distancia crítica en kiloparsecs:
# r = sqrt(G * M^2 * t^2 / (v * M)) = sqrt(G * M / v) en metros
#distancia_critica_metros = math.sqrt(G * masa_grande / umbral_velocidad) #
distancia_critica_kpc = 0.2 # 652,313 años luz

global exc
exc=0
def calcular_fuerza(masa1, masa2, distancia, einstein=False):
    """
    Calcula la fuerza gravitatoria entre dos estrellas utilizando la ley de gravitación de Newton o
    una aproximación simplificada basada en la relatividad general.

    Parámetros:
        masa1 (float): Masas solares.
        masa2 (float): Masas solares.
        distancia (float): En kpc.
        einstein (bool, opcional): Si es True, utiliza la aproximación basada en la relatividad general.

    Retorna:
        float: Fuerza gravitatoria entre las dos estrellas en newtons.
    """
    # Constante de gravitación universal (m^3/kg/s^2)
    constante_gravitacion = 6.67430e-11
    # Velocidad de la luz (m/s)
    velocidad_luz = 299792458.0
    
    #convertir unidades
    masa1_kg = masa1 * 1.9885e30  # Convertir masas solares a kg
    masa2_kg = masa2 * 1.9885e30
    distancia_metros = distancia * 3.086e19  # kpc a metros
    
    global exc
    
    if einstein:   
        # Fórmula basada en la relatividad general para la "fuerza equivalente"
        fuerza = (masa1_kg * masa2_kg * constante_gravitacion) / (distancia_metros**2 * velocidad_luz**2)
    else:
        # Fórmula clásica de la ley de gravitación de Newton
        fuerza = 0.0
        try:
            fuerza = constante_gravitacion * (1 / distancia_metros**2) * masa2_kg * masa1_kg
        except:
            exc += 1
            print(f"excepcion [{exc}] m1:{masa1:.2e} m2:{masa2:.2e}")
        # fuerza en Newtons

    return fuerza


def velocidad(grafo):
    estrellas = grafo.nodes
    num_estrellas = len(estrellas)

    # Constantes
    masa_del_sol_kg = 1.9885e30  # Masa del Sol en kilogramos
    kpc_a_metros = 3.086e19  # Kiloparsec a metros

    # Iterar a través de las estrellas
    contador = 0
    porcentaje_completado = 0
        
    for i, estrella in enumerate(estrellas):
        for j, estrella_grande in enumerate(estrellas):
            if j >= i:  # Solo calcular para las estrellas por encima de la diagonal
                continue
                    # Calcula la distancia entre las estrellas estrella1 y estrella2
            distancia_kpc = np.sqrt \
                        (   (estrellas[estrella]['x'] - estrellas[estrella_grande]['x'])**2 + 
                            (estrellas[estrella]['y'] - estrellas[estrella_grande]['y'])**2 + 
                            (estrellas[estrella]['z'] - estrellas[estrella_grande]['z'])**2)

            # Comprobar si la distancia cumple con el criterio
            if distancia_kpc > distancia_critica_kpc:
                continue
            
            estrellas[estrella]['entorno'] += 1
            estrellas[estrella_grande]['entorno'] += 1
                        
            Masa_sol = masa_del_sol_kg
            M1 = estrellas[estrella]['masa']
            M2 = estrellas[estrella_grande]['masa']
            c_kpc_a_metros = kpc_a_metros
            fuerza = calcular_fuerza(M1, M2, distancia_kpc)  # Fuerza en newtons
            m1_x = estrellas[estrella]['x']*kpc_a_metros
            m1_y = estrellas[estrella]['y']*kpc_a_metros
            m1_z = estrellas[estrella]['z']*kpc_a_metros
            m2_x = estrellas[estrella_grande]['x']*kpc_a_metros
            m2_y = estrellas[estrella_grande]['y']*kpc_a_metros
            m2_z = estrellas[estrella_grande]['z']*kpc_a_metros
            
            # Calcular componentes de la fuerza en x, y, z
            fuerza_x = fuerza * abs(m1_x - m2_x) / (c_kpc_a_metros * distancia_kpc) 
            fuerza_y = fuerza * abs(m1_y - m2_y) / (c_kpc_a_metros * distancia_kpc) 
            fuerza_z = fuerza * abs(m1_z - m2_z) / (c_kpc_a_metros * distancia_kpc)

            # La velocidad que aporta es la fuerza dividida por su masa
            # y pasamos estas velocidades de m/seg a km/seg
            # Actualizar velocidades de las estrellas
            estrellas[estrella]['vx_t1'] += (fuerza_x / (M1*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s
            estrellas[estrella]['vy_t1'] += (fuerza_y / (M1*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s
            estrellas[estrella]['vz_t1'] += (fuerza_z / (M1*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s

            # Aplicar la acción inversa en la estrella grande
            estrellas[estrella_grande]['vx_t1'] -= (fuerza_x / (M2*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s
            estrellas[estrella_grande]['vy_t1'] -= (fuerza_y / (M2*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s
            estrellas[estrella_grande]['vz_t1'] -= (fuerza_z / (M2*masa_del_sol_kg)) / 1000.0  # Dividido por 1000 para convertir a km/s            
        
        
        
        contador += 1  # Incrementar el contador

        # Calcular el porcentaje completado
        test = porcentaje_completado
        porcentaje_completado = ((contador / num_estrellas) * 100) // 1

        # Imprimir el porcentaje en la misma línea
        if test != porcentaje_completado: 
            print(f"Progreso: {porcentaje_completado:.2f}%, ultima entorno:{estrellas[estrella]['entorno']}")

        # Esperar un momento para simular el progreso
        # time.sleep(0.1)  # Puedes ajustar este valor para controlar la velocidad de actualización

    # Calcular aceleraciones y actualizar velocidades para cada estrella
    for i, estrella in enumerate(estrellas):
        # Diferencia entre la velocidad teórica y la velocidad efectiva (GSM)
        estrellas[estrella]['dif_vx'] = estrellas[estrella]['vx'] - estrellas[estrella]['vx_t1']
        estrellas[estrella]['dif_vy'] = estrellas[estrella]['vy'] - estrellas[estrella]['vy_t1']
        estrellas[estrella]['dif_vz'] = estrellas[estrella]['vz'] - estrellas[estrella]['vz_t1']
        
    guardar_grafo_a_csv(grafo, grafo_csv)

    # Imprimir mensaje de finalización
    print("\nCálculo completado.")




archivo_grafo = 'datos/grafo_estrellas.gpickle'
archivo_copia = 'datos/grafo_estrellas_backup.gpickle'
grafo_csv = 'datos/grafo.csv'

# Comprobar si el archivo del grafo existe
try:
    with open(archivo_grafo, 'rb') as f:
        grafo = pickle.load(f)
        print(f"1 - procesado binario {archivo_grafo}, {grafo.number_of_nodes()} estrellas")
        #comprobar si faltan calculos y completarlos
        grafo=evolucionar_grafo(grafo, True)     
        guardar_grafo_a_csv(grafo, grafo_csv)
                               
    # Crear una copia de seguridad del archivo del grafo
    shutil.copyfile(archivo_grafo, archivo_copia)
     
except FileNotFoundError:
    
    try:
        grafo=cargar_grafo_desde_csv(grafo_csv)
        print(f"2 - procesado grafo_csv {grafo_csv}, {grafo.number_of_nodes()} estrellas")        
        
    except FileNotFoundError:
        # si nada es posible entoces crea un grafo nuevo a partir de los datos de Gaia
        print("El archivo del grafo no existe. Creando el grafo...")
        grafo = crear_grafo()
        velocidad(grafo)
        print(f"3 - creado nuevo, {grafo.number_of_nodes()} estrellas")                

actualizar_grafo_en_disco(grafo, archivo_grafo)

ver(grafo)
print("Trabajo terminado")
