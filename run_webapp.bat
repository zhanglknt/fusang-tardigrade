@echo off
REM ============================================================
REM Fusang Web Server Launcher (Windows)
REM ============================================================
REM Automatically activates conda env and starts the web server.
REM ============================================================

setlocal enabledelayedexpansion

set ENV_NAME=fusang-web
set SCRIPT_DIR=%~dp0

REM --- Activate conda and run ---
call conda activate %ENV_NAME% 2>nul || (
    echo [ERROR] Cannot activate conda environment '%ENV_NAME%'.
    echo Run 'setup_conda.bat' first to create the environment.
    exit /b 1
)

echo.
echo ============================================================
echo   Fusang Web Server
echo   Access at: http://localhost:5000
echo ============================================================
echo.

cd /d "%SCRIPT_DIR%"
python fusang_webapp.py %*
