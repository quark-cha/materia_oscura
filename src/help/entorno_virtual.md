# Cómo crear y manejar un entorno virtual en Python

## 1. Crear un entorno virtual

✅ Asegúrate de tener Python instalado:

```bash
python3 --version
```

> En Windows puede ser `python` en lugar de `python3`.

✅ Crear el entorno virtual:

```bash
python3 -m venv nombre_del_entorno
```

📦 Esto creará una carpeta llamada `nombre_del_entorno` con todo lo necesario para tu entorno.

Ejemplo:

```bash
python3 -m venv venv
```

---

## 2. Activar el entorno virtual

🔵 En Linux / WSL / macOS:

```bash
source venv/bin/activate
```

🔵 En Windows (CMD):

```cmd
venv\Scripts\activate.bat
```

🔵 En Windows (PowerShell):

```powershell
venv\Scripts\Activate.ps1
```

> Una vez activado, verás algo como `(venv)` al principio de tu línea de comandos.

---

## 3. Instalar paquetes dentro del entorno

Cuando el entorno está activado:

```bash
pip install nombre_paquete
```

Ejemplo:

```bash
pip install matplotlib jupyter pandas
```

> Estos paquetes solo se instalarán en este entorno, no en todo el sistema.

---

## 4. Guardar las dependencias

Después de instalar lo que necesites:

```bash
pip freeze > requirements.txt
```

> Esto genera un archivo con la lista de dependencias.

---

## 5. Restaurar dependencias (en otro entorno o máquina)

```bash
pip install -r requirements.txt
```

---

## 6. Salir del entorno virtual

Cuando termines:

```bash
deactivate
```

---

## 7. Eliminar el entorno virtual (opcional)

Solo borra la carpeta:

```bash
rm -rf venv
```

> o en Windows, elimínala desde el Explorador de archivos.
