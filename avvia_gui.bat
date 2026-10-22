@echo off
title Energyinator Studio
echo ===================================================
echo     Avvio di Energyinator Studio (Streamlit GUI)...
echo ===================================================
echo.
cd /d "%~dp0"
.\.venv\python.exe -m streamlit run app_gui.py
pause

