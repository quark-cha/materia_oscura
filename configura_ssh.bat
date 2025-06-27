@echo off
setlocal

echo ===============================
echo Verificando clave SSH existente...
echo ===============================
set KEY=%USERPROFILE%\.ssh\id_ed25519

IF EXIST "%KEY%" (
    echo Ya existe una clave SSH: %KEY%
) ELSE (
    echo Generando nueva clave SSH...
    ssh-keygen -t ed25519 -C "tu_correo@ejemplo.com" -f "%KEY%" -N ""
)

echo.
echo Copiando tu clave pública al portapapeles...
type "%KEY%.pub" | clip

echo.
echo ✅ Clave pública copiada al portapapeles.
echo 💡 Ahora ve a GitHub y pégala aquí:
start https://github.com/settings/keys

pause

echo.
echo Probando conexión SSH con GitHub...
ssh -T git@github.com

pause
endlocal
