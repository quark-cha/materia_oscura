# src/dependencias.py
"""
Verificación esencial de dependencias científicas
"""

# Lista de dependencias críticas (como tuplas)
DEPENDENCIAS = [
    ('pandas', 'pandas==2.0.3'),
    ('numpy', 'numpy==1.24.4'),
    ('matplotlib', 'matplotlib==3.7.2'),
    ('networkx', 'networkx==3.1'),
    ('astropy', 'astropy==5.3'),
    ('scipy', 'scipy==1.10.1'),
    ('pyvo', 'pyvo==1.4'),
    ('seaborn', 'seaborn==0.13.2'),
    ('pathlib', 'pathlib')
]

def verificar_dependencias():
    """Verifica si todas las dependencias están instaladas"""
    import sys
    import importlib.util
    
    problemas = False
    
    # CORRECCIÓN: Usar índices de tupla en lugar de claves de diccionario
    for dep in DEPENDENCIAS:
        mod_name = dep[0]  # Primer elemento: nombre del módulo
        pip_spec = dep[1]  # Segundo elemento: comando pip
        
        # Verificar si el módulo existe
        if not importlib.util.find_spec(mod_name):
            print(f"❌ ERROR: Módulo científico faltante: {mod_name}")
            print(f"    Solución: pip install {pip_spec}")
            problemas = True
    
    if problemas:
        print("\n" + "="*70)
        print("POR FAVOR INSTALE LAS DEPENDENCIAS FALTANTES Y EJECUTE NUEVAMENTE")
        print("="*70)
        sys.exit(1)