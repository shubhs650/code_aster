@echo off
setlocal EnableExtensions

REM Wrapper for the packaged Windows layout.
REM Installed next to run_aster.bat in code_aster's bin directory.

chcp 65001 >NUL
set "PYTHONIOENCODING=utf-8"
set "LC_ALL=en_US.UTF-8"
set "LANG=en_US.UTF-8"
set "LANGUAGE=en_US.UTF-8"

for %%# in ("%~dp0..") do set "RUNASTER_ROOT=%%~f#"
for %%# in ("%RUNASTER_ROOT%\..") do set "PACKAGED_DIR=%%~f#"

set "PYTHON_DIR=%PACKAGED_DIR%\Python3.10"
set "EXTERNAL_DIR=%PACKAGED_DIR%\External"
set "CMAKE_DIR=%EXTERNAL_DIR%\CMake\bin"

set "ASTER_LIBDIR=%RUNASTER_ROOT%\lib\aster"
set "ASTER_DATADIR=%RUNASTER_ROOT%\share\aster"
set "ASTER_LOCALEDIR=%RUNASTER_ROOT%\share\locale\aster"
set "ASTER_ELEMENTSDIR=%RUNASTER_ROOT%\lib\aster"

set "MEDCOUPLING_PY=%EXTERNAL_DIR%\medcoupling-9.11.0\lib\python3.10\site-packages"

echo [run_aster] RUNASTER_ROOT=%RUNASTER_ROOT%
echo [run_aster] PACKAGED_DIR=%PACKAGED_DIR%
echo [run_aster] PYTHON_DIR=%PYTHON_DIR%
echo [run_aster] EXTERNAL_DIR=%EXTERNAL_DIR%
echo [run_aster] ASTER_LIBDIR=%ASTER_LIBDIR%
echo [run_aster] ASTER_DATADIR=%ASTER_DATADIR%
echo [run_aster] MEDCOUPLING_PY=%MEDCOUPLING_PY%

set "PYTHONPATH=%ASTER_LIBDIR%;%MEDCOUPLING_PY%;%PYTHONPATH%"
set "LD_LIBRARY_PATH=%PYTHON_DIR%;%EXTERNAL_DIR%\hdf51.10.5\bin;%EXTERNAL_DIR%\hdf51.10.5\lib;%EXTERNAL_DIR%\MED-4.4.1\lib;%EXTERNAL_DIR%\mfront-4.2.3\bin;%EXTERNAL_DIR%\mgis-2.2\bin;%EXTERNAL_DIR%\medcoupling-9.11.0\lib;%MEDCOUPLING_PY%;%EXTERNAL_DIR%;%ASTER_LIBDIR%"
set "PATH=%RUNASTER_ROOT%\bin;%PYTHON_DIR%;%ASTER_LIBDIR%;%EXTERNAL_DIR%;%CMAKE_DIR%;%EXTERNAL_DIR%\hdf51.10.5\bin;%EXTERNAL_DIR%\hdf51.10.5\lib;%EXTERNAL_DIR%\MED-4.4.1\lib;%EXTERNAL_DIR%\mfront-4.2.3\bin;%EXTERNAL_DIR%\mgis-2.2\bin;%EXTERNAL_DIR%\medcoupling-9.11.0\lib;%MEDCOUPLING_PY%;%PATH%"

echo [run_aster] PYTHONPATH=%PYTHONPATH%
echo [run_aster] PATH=%PATH%

if not exist "%PYTHON_DIR%\python.exe" (
    echo [ERROR] Missing Python interpreter: "%PYTHON_DIR%\python.exe"
    endlocal & exit /b 1
)

call "%PYTHON_DIR%\python.exe" -m run_aster.run_aster_main %*
set "RET=%ERRORLEVEL%"

endlocal & exit /b %RET%
