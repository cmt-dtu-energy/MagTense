@echo off
REM ---------------------------------------------------------------------------
REM  MagTense Windows installer launcher.
REM
REM  Double-click this file (or run it from a terminal) to install a complete
REM  MagTense development environment. It simply runs install-magtense.ps1 with
REM  an execution policy that allows the script to run for this process only,
REM  so you do not have to change any system-wide PowerShell settings.
REM
REM  Any arguments you pass are forwarded to the PowerShell script, e.g.:
REM      install-magtense.bat -Compute cpu -PyVersion 313
REM ---------------------------------------------------------------------------

setlocal
set "PS1=%~dp0install-magtense.ps1"

if not exist "%PS1%" (
    echo Could not find install-magtense.ps1 next to this launcher.
    echo Make sure both files are in the same folder.
    pause
    exit /b 1
)

powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%PS1%" %*
set "RC=%ERRORLEVEL%"

echo.
if not "%RC%"=="0" (
    echo Installation reported an error ^(exit code %RC%^). See the log path above.
)
echo Press any key to close this window.
pause >nul
exit /b %RC%
