import sys
import subprocess
import platform

def verificar_entorno():
    print("🔍 Iniciando diagnóstico científico...")
    print(f"Sistema: {platform.system()} {platform.release()}")
    print(f"Python: {sys.version}")
    
    try:
        import pip
        print("✓ pip instalado")
    except ImportError:
        print("❌ pip no disponible - Instalando...")
        subprocess.check_call([sys.executable, "-m", "ensurepip", "--default-pip"])
    
    dependencias = ["pandas", "numpy", "matplotlib", "astropy"]
    problemas = []
    
    for paquete in dependencias:
        try:
            mod = __import__(paquete)
            print(f"✓ {paquete} ({mod.__version__}) instalado")
        except ImportError:
            print(f"❌ {paquete} faltante")
            problemas.append(paquete)
    
    if problemas:
        print("\n🚨 PROBLEMAS DETECTADOS - INICIANDO REPARACIÓN...")
        subprocess.check_call([sys.executable, "-m", "pip", "install"] + problemas)
        print("✅ Reparación completada - Por favor vuelve a ejecutar tu análisis")
    else:
        print("\n✅ Entorno científico VERIFICADO - Todo funciona correctamente")
        print("Puedes ejecutar: python src\grafo-galaxia.py")

if __name__ == "__main__":
    verificar_entorno()