@echo off
REM build-sidecar.bat — builds the Mathilda sidecar binary for Windows.
REM Requires MSYS2 MinGW64 (gcc, make, mingw-w64-x86_64-gmp) in PATH.
REM Run this script from a standard PowerShell or CMD prompt; MSYS2 tools
REM must already be on PATH (e.g. C:\msys64\mingw64\bin and C:\msys64\usr\bin).
setlocal enabledelayedexpansion

set "SCRIPT_DIR=%~dp0"
REM Strip trailing backslash
if "%SCRIPT_DIR:~-1%"=="\" set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"
set "REPO_ROOT=%SCRIPT_DIR%\.."
set "BINARIES_DIR=%SCRIPT_DIR%\src-tauri\binaries"

REM Determine the Rust target triple for naming the sidecar binary.
for /f "tokens=2 delims= " %%i in ('rustc -vV 2^>nul ^| findstr /C:"host:"') do set "TARGET=%%i"
if "%TARGET%"=="" (
    echo ERROR: Could not determine Rust target triple.
    echo        Ensure rustc is installed and on PATH.
    exit /b 1
)

echo Building Mathilda for target: %TARGET%

REM Build using make. GMP is provided by MSYS2 mingw-w64-x86_64-gmp.
REM Readline is not available on Windows, so we pass USE_READLINE=0.
REM ECM vendored build is skipped (USE_ECM=0) to keep the Windows build simple.
pushd "%REPO_ROOT%"
make USE_ECM=0 USE_READLINE=0 USE_GRAPHICS=0 USE_LAPACK=0 -j4
if errorlevel 1 (
    echo ERROR: make failed. Check that gcc and gmp are in PATH.
    popd
    exit /b 1
)
popd

if not exist "%BINARIES_DIR%" mkdir "%BINARIES_DIR%"

set "SRC_EXE=%REPO_ROOT%\Mathilda.exe"
set "DST_EXE=%BINARIES_DIR%\mathilda-%TARGET%.exe"

if not exist "%SRC_EXE%" (
    echo ERROR: Build succeeded but Mathilda.exe not found at %SRC_EXE%
    exit /b 1
)

copy /y "%SRC_EXE%" "%DST_EXE%"
if errorlevel 1 (
    echo ERROR: Failed to copy sidecar to %DST_EXE%
    exit /b 1
)

echo.
echo Sidecar installed: %DST_EXE%
echo.
echo Next steps:
echo   cd %SCRIPT_DIR%
echo   cargo tauri build
endlocal
