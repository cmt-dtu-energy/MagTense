@ECHO OFF

setlocal
pushd %~dp0

REM Command file for the MagTense TechManual.
REM
REM   make html    build the manual and publish it into ..\ (the docs folder
REM                that GitHub Pages serves), then commit the result
REM   make check   build only, leaving the published site untouched
REM
REM The publish step runs only after a successful build. An earlier version of
REM this file deleted ..\_images, ..\_sources and ..\_static before checking
REM whether Sphinx had produced anything, so a failed build silently wiped the
REM published site.

REM Invoke Sphinx through the active interpreter rather than through whichever
REM sphinx-build.exe happens to come first on PATH. conf.py imports
REM sphinx_rtd_theme, so picking up a different environment's sphinx-build
REM fails the build with "No module named 'sphinx_rtd_theme'".
if "%SPHINXBUILD%" == "" (
	set "SPHINXBUILD=python -m sphinx"
)
set "SOURCEDIR=."
set "BUILDDIR=_build"
set "PUBLISHDIR=.."

if "%1" == "" goto help
if /i "%1" == "help" goto help
if /i "%1" == "check" goto check
if /i "%1" == "preview" goto preview

%SPHINXBUILD% --version >NUL 2>NUL
if errorlevel 1 goto nosphinx

%SPHINXBUILD% -M %1 "%SOURCEDIR%" "%BUILDDIR%" %SPHINXOPTS% %O%
if errorlevel 1 goto buildfailed

REM Only the html target produces something that can be published.
if /i not "%1" == "html" goto end

if not exist "%BUILDDIR%\html\index.html" (
	echo.
	echo.No index.html was produced - %PUBLISHDIR% left untouched.
	goto fail
)

echo.
echo.Publishing into %PUBLISHDIR% ...

REM Delete the previously published output first, so that pages whose source
REM has been removed do not linger on the website. README.md and .nojekyll are
REM hand-maintained and are deliberately not touched.
if exist "%PUBLISHDIR%\_images"  rmdir /s /q "%PUBLISHDIR%\_images"
if exist "%PUBLISHDIR%\_sources" rmdir /s /q "%PUBLISHDIR%\_sources"
if exist "%PUBLISHDIR%\_static"  rmdir /s /q "%PUBLISHDIR%\_static"
del /q "%PUBLISHDIR%\*.html" >NUL 2>NUL
del /q "%PUBLISHDIR%\*.js" >NUL 2>NUL
del /q "%PUBLISHDIR%\*.inv" >NUL 2>NUL
del /q "%PUBLISHDIR%\.buildinfo" >NUL 2>NUL

xcopy "%BUILDDIR%\html" "%PUBLISHDIR%" /e /i /y /q >NUL
if errorlevel 1 (
	echo.
	echo.Copying the built pages into %PUBLISHDIR% failed.
	goto fail
)

rmdir /s /q "%BUILDDIR%"

echo.
echo.Done. Commit the changes under docs\ to update the website.
goto end

:preview
REM Build for inspection only. Nothing outside _build is touched, so this is
REM always safe to run, and warnings do not stop it producing output.
%SPHINXBUILD% --version >NUL 2>NUL
if errorlevel 1 goto nosphinx
%SPHINXBUILD% -b html "%SOURCEDIR%" "%BUILDDIR%\html" %SPHINXOPTS% %O%
if errorlevel 1 goto buildfailed
echo.
echo.Open this file in a browser:
echo.
echo.    %~dp0_build\html\index.html
echo.
goto end

:check
%SPHINXBUILD% --version >NUL 2>NUL
if errorlevel 1 goto nosphinx
%SPHINXBUILD% -b html -W -n "%SOURCEDIR%" "%BUILDDIR%\html" %SPHINXOPTS% %O%
if errorlevel 1 goto buildfailed
echo.
echo.Build is clean. Nothing was published - run "make html" for that.
goto end

:nosphinx
echo.
echo.Sphinx was not found in the active environment.
echo.
echo.Install it with:
echo.
echo.    pip install -r ..\requirements.txt
echo.
echo.or set the SPHINXBUILD environment variable to the sphinx-build
echo.executable you want to use.
echo.
goto fail

:buildfailed
echo.
echo.The Sphinx build FAILED - nothing was published and %PUBLISHDIR% is unchanged.
echo.
goto fail

:help
%SPHINXBUILD% -M help "%SOURCEDIR%" "%BUILDDIR%" %SPHINXOPTS% %O%
echo.
echo.MagTense targets:
echo.  preview  build into _build\html for inspection, publishing nothing
echo.  check    build with warnings as errors, publishing nothing
echo.  html     build and publish into the docs folder
goto end

:fail
popd
endlocal
exit /b 1

:end
popd
endlocal
