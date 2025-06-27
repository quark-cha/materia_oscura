## 🚀 Configuración Rápida

Este proyecto está diseñado para ejecutarse directamente con las dependencias del sistema. **No se requiere entorno virtual** para simplificar la reproducción.
No obstante puedes querer tener uno para que no interfiera con otros proyectos tuyos, en ese caso, si no lo tienes has de crear el manejo de este entorno:
(Opcionalmente)
Para ver como hacerlo lee el documento (Materia_oscura)src/help/entorno_virtual.md


### Requisitos previos
- Python 3.8+
- PIP (gestor de paquetes)

### Instalación en 2 pasos:
1. Situate en el directorio donde quieras tener este u otros trabajos
   tus proyectos repuestos colgaran de ahí cuando clones los repositorios
2. Clona el repositorio:
   ```bash
   git clone https://github.com/quark-cha/materia_oscura.git
   cd materia_oscura
3. Si sólo quieres esplorar el trabajo sin instalar nada usa el navegador
   y abre el link: https://github.com/quark-cha/materia_oscura


Sistema de Análisis Astrofísico - Victor Estrada Díaz

Lugar oficial de todos mis trabajos y teorias: estradad.es/teorias
Zenodo: https://zenodo.org/uploads/15742442
GitHub: quark-cha/materia_oscura

---

Estructura del Directorio
-------------------------

src/
├── versiones/               # Zona sagrada de trabajos históricos
│   ├── v3/                  # Versión exacta usada en Zenodo 15733011
│   │   
│   └── v4/                  # Nueva versión (evolución)
│       └── grafo-galaxia.py # Bloques 1-3 originales + 11-13 nuevos
│
├── grafo-galaxia.py         # Versión activa (actualmente v4)
├── test_reproducibilidad.py # Verificador científico (pendiente)
└── memoria_evolucion.md     # Registro de cambios y justificaciones

Directorios de trabajo:
- Figuras: /plots/tmp
- Datos: /datos
- documentos: /doc   # Resultados las investigaciones y desarroyo de ideas 

---

Cómo Usar el Sistema
--------------------

1. Reproducir resultados históricos (ej. v3)

Ejecutar versión exacta de un trabajo publicado:
python versiones/v3/grafo-galaxia.py

Verificar reproducibilidad científica:
python test_reproducibilidad.py --version 3

2. Ejecutar versión actual con mejoras

Usar nueva visualización (bloque 12) manteniendo cálculos históricos:
python grafo-galaxia.py --version 4

3. Añadir nuevas funcionalidades
- Crear nuevo bloque numérico (n > último existente)
- Implementar sin modificar bloques anteriores
- Registrar cambios en memoria_evolucion.md
- Verificar compatibilidad:
  python test_reproducibilidad.py --all

---

Bloques Clave Actuales

Históricos (inmutables)
Bloque  Función                 Versión
1       Carga datos Gaia        v3
2       Construcción grafo base v3
3       Cálculo propiedades     v3

Nuevos (evolución)
Bloque  Función                     Versión
11      Verificador reproducibilidad v4
12      Visualización mejorada       v4
13      Sistema ejecución versionada v4

---

Flujo de Trabajo
----------------

1. Para análisis histórico:
   - Ejecutar: /src/versiones/v3/grafo-galaxia.py
   - Resultados en /plots/tmp (figuras) y /datos (csv)

2. Para desarrollo nuevo:
   - Añadir bloques con números superiores
   - Probar con test_reproducibilidad.py
   - Documentar en memoria_evolucion.md

3. Para verificación científica:
   python test_reproducibilidad.py --version 3
   python test_reproducibilidad.py --all

---

Documentación Relacionada
-------------------------
Explicación teórica: estradad.es/teorias
Base de datos Gaia: cosmos.esa.int/web/gaia
Resultados clave en Zenodo: zenodo.org/records/15733011

---

"La ciencia avanza cuando los datos obligan a cambiar las teorías,  
no cuando las teorías obligan a ignorar los datos."  
- Victor Estrada Díaz