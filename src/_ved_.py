
import sys
from pathlib import Path

# AÑADIR LA RAÍZ DEL PROYECTO A sys.path
ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))
    
import matplotlib
import matplotlib.pyplot as plt
import math

import networkx as nx
import os
import numpy as np
import seaborn as sns
from scipy.interpolate import CubicSpline
from mpl_toolkits.axes_grid1 import make_axes_locatable

global tipo
tipo = ""

print(matplotlib.__version__)

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

def milist(var, var_name, n=5):
    print(f" {var_name} (primeros {n} elementos): {{")
    for i in range(0,n):
        print(f"  '{i}': {var[i]},")
    print("}")
    

    
def ver(G):

    #PLOT función para poner un grafico y asignar un numero de figura

    # Definir las variables para la posición y velocidad del Sol en el sistema de coordenadas galácticas
    x0 = 8.3 # kpc
    y0 = 0 # kpc
    z0 = 0.027 # kpc
    vx0 = -11.1 # km/s
    vy0 = 232.24 # km/s
    vz0 = 7.25 # km/s
    
    print("Visualizando grafo, tamaño del grafo:", G.number_of_nodes())
    print("atributos grafo", G.node_dict_factory.values)

    # Calcular la velocidad tangencial del Sol en el plano galáctico
    v0 = math.sqrt(vx0**2 + vy0**2)
    d0 = math.sqrt(x0**2 + z0**2)
    

    """
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
            'pmdec': pmdec,         #movimiento aparente en declinación
            'vx': 0.0,              #velocidad real (GSM) en coordenada X
            'vy': 0.0,              #velocidad real (GSM) en coordenada Y
            'vz': 0.0,              #velocidad real (GSM) en coordenada Z
            'entorno': 0,           #numero de estrellas en el entorno de influencia de la estrella
            'dif_vx':0.0,           #diferencia velocidad real - velocidad teorica t1 coordenada X
            'dif_vy':0.0,           #diferencia velocidad real - velocidad teorica t1 coordenada Y
            'dif_vz':0.0,           #diferencia velocidad real - velocidad teorica t1 coordenada Z
            'vx_t1': 0.0,           #velocidad teorica t1 en coordenada X
            'vy_t1': 0.0,           #velocidad teorica t1 en coordenada Y
            'vz_t1': 0.0,           #velocidad teorica t1 en coordenada Z
            'mod_v': 0.0,           #modulo de la velocidad teorica t1
            'mod_t': 0.0,           #modulo de la velocidad tangencial teorica t1
            'divergencia':[0,0,0]   #divergencia velocidad teorica t1 - velocidad real
            
    """
    # Convertir los atributos de los nodos a arrays

    masa = np.array(list(nx.get_node_attributes(G, 'masa').values()))
    
    radios = np.array(list(nx.get_node_attributes(G, 'radio').values()))
    mod_t = np.array(list(nx.get_node_attributes(G, 'mod_t').values()))
    mod_v = np.array(list(nx.get_node_attributes(G, 'mod_v').values()))
    diver = np.array(list(nx.get_node_attributes(G, 'divergencia').values()))

    x = np.array(list(nx.get_node_attributes(G, 'x').values()))
    y = np.array(list(nx.get_node_attributes(G, 'y').values()))
    z = np.array(list(nx.get_node_attributes(G, 'z').values()))
    
    vx = np.array(list(nx.get_node_attributes(G, 'vx').values()))
    vy = np.array(list(nx.get_node_attributes(G, 'vy').values()))
    vz = np.array(list(nx.get_node_attributes(G, 'vz').values()))
    
    vx_t1 = np.array(list(nx.get_node_attributes(G, 'vx_t1').values()))
    vy_t1 = np.array(list(nx.get_node_attributes(G, 'vy_t1').values()))
    vz_t1 = np.array(list(nx.get_node_attributes(G, 'vz_t1').values()))
    
    dif_vx = np.array(list(nx.get_node_attributes(G, 'dif_vx').values()))
    dif_vy = np.array(list(nx.get_node_attributes(G, 'dif_vy').values()))
    dif_vz = np.array(list(nx.get_node_attributes(G, 'dif_vz').values()))
    
    vx_t1 = np.array(list(nx.get_node_attributes(G, 'vx_t1').values()))
    vy_t1 = np.array(list(nx.get_node_attributes(G, 'vy_t1').values()))
    vz_t1 = np.array(list(nx.get_node_attributes(G, 'vz_t1').values()))
    
    entorno = np.array(list(nx.get_node_attributes(G, 'entorno').values()))
    
    
    hecho = []
    figuras = [-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11]
    
    figura = - 1
    print ("Figura ",figura)
    if not figura in hecho and figura in figuras:
        # Crear figura
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        rango_xyz = (-100, 100) #Kpc 
        ax.set_xlim(rango_xyz)
        ax.set_ylim(rango_xyz)
        ax.set_zlim(rango_xyz)
        
        # Título y leyenda
        ax.set_title('Situacion de la muestra del grafo')

        # Gráfico 3D de las coordenadas x, y, z
        ax.scatter(x, y, z, color='black', s=0.5)

        # Marcar un punto del sol
        ax.scatter(x0, y0, z0, c='r', marker='o', s=10, zorder=10)

        plot(fig,figura)

#----------------------------------------------------------------
    figura -=1
    print ("Figura ",figura)
    if not figura in hecho and figura in figuras:

        # Crear figura
        fig, ax = plt.subplots()  
      
        # Scatter plot
        sc=ax.scatter(radios, mod_t, s=0.1, zorder=9, c='green')
      
        # Título y leyenda
        ax.set_title('Modulo velocidad tangencial real GSM')
        
        # Agregar un subtitulo a la derecha utilizando fig.text()
        fig.text(0.85, 0.5, 'v.real=v.aparente-v.sol', va='center', rotation=90, fontsize=8)
         
        # Marcar un punto del sol
        ax.scatter([d0], [v0], c='r', marker='o', s=10, zorder=10)
    
        ax.set_xlabel('Distancia al Centro Galáctico (kpc)')
        ax.set_ylabel('Mod_t (km/s)')

        # Establecer límites de los ejes
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 1200)
    
        # Mostrar la figura
        plot(fig,figura)

#----------------------------------------------------------------
    figura -=1
    print ("Figura ",figura)
    if not figura in hecho and figura in figuras:

        # Crear figura
        fig, ax = plt.subplots()  
      
        # Scatter plot
        sc=ax.scatter(radios, mod_v, s=0.1, zorder=9, c='blue')
      
        # Título y leyenda
        ax.set_title('Modulo velocidad real GSM')
        # Agregar un subtitulo a la derecha utilizando fig.text()
        fig.text(0.85, 0.5, 'v.real=v.aparente-v.sol', va='center', rotation=90, fontsize=8)
        
        # Marcar un punto del sol
        ax.scatter([d0], [v0], c='r', marker='o', s=10, zorder=10)
    
        ax.set_xlabel('Distancia al Centro Galáctico (kpc)')
        ax.set_ylabel('Mod_v (km/s)')

        # Establecer límites de los ejes
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 1200)
    
        # Mostrar la figura
        plot(fig,figura)

#----------------------------------------------------------------

    figura -=1
    print ("Figura ",figura)
    
    if not figura in hecho and figura in figuras:

        # Crear figura
        fig, ax = plt.subplots()  
      
        # Scatter plot
        sc=ax.scatter(radios, diver, s=0.1, zorder=9, c='blue',label='divergencia')
      
        # Título y leyenda
        ax.set_title('Modulo de divergencia velocidad real GSM')
        # Agregar un subtitulo a la derecha utilizando fig.text()
        fig.text(0.85, 0.5, 'v.divergente=v.real-v.tangencial', va='center', rotation=90, fontsize=8)
        
        # Marcar un punto del sol
        ax.scatter([d0], [v0], c='r', marker='o', s=10, zorder=10)
    
        ax.set_xlabel('Distancia al Centro Galáctico (kpc)')
        ax.set_ylabel('Divergencia (km/s)')

        # Establecer límites de los ejes
        ax.set_xlim(0, 100)
        ax.set_ylim(-140, 1200)
    
        # Mostrar la figura
        plot(fig,figura)

#----------------------------------------------------------------

    figura -=1
    print ("Figura ",figura)
    if not figura in hecho and figura in figuras:

        # Crear figura
        fig, ax = plt.subplots()  

        sc=ax.scatter(masa, mod_v, s=0.1, zorder=2, c='green', label='Masa-velocidad real')

        ax.set_title('Masa relacionada con la velocidad real')

        ax.set_xlabel('Masa en soles')
        ax.set_ylabel('Modulo velocidad real(km/s)')
        
        # Establecer límites de los ejes
        ax.set_xlim(0, 1000)
        ax.set_ylim(-10, 2500)

        ax.legend()
        
        # Mostrar la figura
        plot(fig,figura)

#----------------------------------------------------------------

    figura -=1
    print ("Figura ",figura)
    if not figura in hecho and figura in figuras:

        
        # Crear figura
        fig, ax = plt.subplots()

        sc=ax.scatter(masa, mod_t, s=0.1, zorder=2, c='green', label='Masa-velocidad tangencial')

        ax.set_title('Masa relacionada con la velocidad tangencial')

        ax.set_xlabel('Masa en soles')
        ax.set_ylabel('Modulo velocidad tangencial(km/s)')
        ax.legend()
        
        # Establecer límites de los ejes
        ax.set_xlim(0, 1000)
        ax.set_ylim(-10, 2500)


        # Mostrar la figura
        plot(fig,figura)

#----------------------------------------------------------------

    figura -=1
    print ("Figura ",figura)
    if not figura in hecho and figura in figuras:

        
        # Crear figura
        fig, ax = plt.subplots()
        
        ax.set_title('Entorno de influencia de las estrellas')

        sc=ax.scatter(radios, entorno, c=masa, 
                      vmax=200.0, cmap='plasma', s=0.1, zorder=2, label='Entorno')
        # Agregar barra de color debajo del eje y
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(sc, cax=cax, label='Masa de las estrellas (masas solares)')


        ax.set_xlabel('Distancia al centro galactico en Kpc')
        ax.set_ylabel('estrellas cercanas influyentes')
        ax.legend()
        
        # Establecer límites de los ejes
        ax.set_xlim(0, 8)
        ax.set_ylim(0, 3000)


        # Mostrar la figura
        plot(fig,figura)

#----------------------------------------------------------------

    figura -=1
    print ("Figura ",figura)
    if not figura in hecho and figura in figuras:

        
        # Crear figura
        fig, ax = plt.subplots()

        ax.set_title('Velocidad teorica t1')

        milist(vx_t1,"vx_t1",5)

        sc=ax.scatter(radios, vx_t1, s=0.1, zorder=2, c='green', label='vxt1')
        sc=ax.scatter(radios, vy_t1, s=0.1, zorder=2, c='red', label='vyt1')
        sc=ax.scatter(radios, vz_t1, s=0.1, zorder=2, c='blue', label='vzt1')

        # Agregar líneas de color en las etiquetas
        legend = ax.legend()
        for handle in legend.legend_handles:
            handle.set_sizes([30])  # Ajusta el tamaño de la línea de color
    

        ax.set_xlabel('Velocidad teorica t1')
        ax.set_ylabel('(km/s)')
        ax.legend()
        
        # Establecer límites de los ejes
        ax.set_xlim(0, 60)
        ax.set_ylim(-2000, 2000)

        # Mostrar la figura
        plot(fig,figura)

#----------------------------------------------------------------


    figura -=1
    print ("Figura ",figura)
    if not figura in hecho and figura in figuras:

        # Crear figura    
        
        fig = plt.figure()
        ax = fig.add_subplot(111)

        rango_x = (0, 60) #Kpc 
        ax.set_xlim(rango_x)

        ax.set_title('Velocidad real x y z')

        sc=ax.scatter(x, vx, s=0.1, zorder=3, c='green',label='Vx g')
        sc=ax.scatter(y, vy, s=0.1, zorder=2, c='blue',label='Vy b')
        sc=ax.scatter(z, vz, s=0.1, zorder=2, c='red',label='Vz r')
        
        # Marcar un punto del sol
        ax.scatter([d0], [v0], c='r', marker='o', s=10, zorder=10)
    
        # Etiquetas personalizadas en la leyenda
        handles, labels = ax.get_legend_handles_labels()
        custom_labels = [f'{label} ({color})' for label, color in zip(labels, ['grey', 'red', 'blue'])]
        ax.legend(handles, custom_labels, loc='upper left')

        ax.set_xlabel('distancia en la galaxia en Kilo prasec')
        ax.set_ylabel('Velocidad en km/s')
        ax.legend()

        # Establecer límites de los ejes
        ax.set_xlim(-60, 60)
        ax.set_ylim(-2000, 2000)

        # Mostrar la figura
        plot(fig,figura)



    figura-=1
    print ("Figura ",figura)
    if not figura in hecho and figura in figuras:

        # Crear figura    
        
        fig = plt.figure()
        ax = fig.add_subplot(111)

        rango_x = (0, 60) #Kpc 
        ax.set_xlim(rango_x)

        ax.scatter(x, vx_t1, s=0.1, zorder=8, c='red')
        ax.scatter(y, vy_t1, s=0.1, zorder=8, c='yellow')
        ax.scatter(z, vz_t1, s=0.1, zorder=8, c='blue')
        # Mostrar la figura
        plot(fig,figura)



#@@##################################################33


def plot(fig, n):
    f = f'{n}'
    titulo = f'Figura {f}'
    plt.suptitle(titulo)
    marca(fig)
    archivo=f'plots/tmp/figura{f}.png'
    try:
        print ("borrando ... ",archivo)
        os.remove(archivo)
    except FileNotFoundError:
        print ("no necesita borrar ... ",archivo)
    fig.savefig(archivo)
    plt.show()
   
    plt.close()
    

def dis(figura,v1,v2,titulo,subtitulo,color,x_etiqueta,y_etiqueta,xlim=0,ylim=0):

    # Crear figura
    fig= plt.figure()
    ax=fig.add_subplot(111)
    
    # Scatter plot
    sc=ax.scatter(v1, v2, s=0.1, zorder=9, c=color)
      
    # Título y leyenda
    ax.set_title(titulo)
        
    # Agregar un subtitulo a la derecha utilizando fig.text()
    fig.text(0.85, 0.5, subtitulo, va='center', rotation=90, fontsize=8)
         
    ax.set_xlabel(x_etiqueta)
    ax.set_ylabel(y_etiqueta)

    # Establecer límites de los ejes
    xmax=np.max(v1)
    ymax=np.max(v2)
    print (" xmax: ", xmax, "ymax: ", ymax)
    
    """
    if not xlim == 0:
        ax.set_xlim(0, xlim)
    else:
        ax.set_xlim=(0,xmax)
    if not ylim == 0:
        ax.set_ylim(0, ylim)
    else:
        ax.set_ylim=(0, ymax)
    """
    # Mostrar la figura
    plot(fig,figura)    