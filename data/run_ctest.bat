@echo off
setlocal EnableExtensions

REM Relocatable CTest launcher for the Packaged Windows layout.
chcp 65001 >NUL
for %%# in ("%~dp0..") do set "RUNASTER_ROOT=%%~f#"
for %%# in ("%RUNASTER_ROOT%\..") do set "PACKAGED_DIR=%%~f#"
set "EXTERNAL_DIR=%PACKAGED_DIR%\External"
set "PYTHON_DIR=%EXTERNAL_DIR%\python-3.10.11"

call "%RUNASTER_ROOT%\share\aster\profile.bat"

if not exist "%EXTERNAL_DIR%\CMake\bin\ctest.exe" (
    echo [ERREUR] CTest absent : "%EXTERNAL_DIR%\CMake\bin\ctest.exe"
    endlocal & exit /b 1
)

"%PYTHON_DIR%\python.exe" -m run_aster.run_ctest_main %*
set "RET=%ERRORLEVEL%"
endlocal & exit /b %RET%
