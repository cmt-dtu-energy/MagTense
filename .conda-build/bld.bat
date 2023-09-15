@echo off
SETLOCAL ENABLEDELAYEDEXPANSION

REM Define dictornary of wheels on PyPI
SET wheels[0].name="2.1.0+cpu_py3.9"
SET wheels[0].link=https://files.pythonhosted.org/packages/56/1a/e7bde7359e07318d610d05ba52484fb5e7318ce693360c1ad6ab26342e22/magtense-2.1.0+cpu-py39-none-win_amd64.whl

SET wheels[1].name="2.1.0+cpu_py3.10"
SET wheels[1].link=https://files.pythonhosted.org/packages/7a/f6/df4567440d519873c780365848e9beca2f59e71ae0818953d872a928c2b5/magtense-2.1.0+cpu-py310-none-win_amd64.whl

SET wheels[2].name="2.1.0+cuda11.6_py3.9"
SET wheels[2].link=https://files.pythonhosted.org/packages/ab/06/7c1b285066e21a0ef19fedb540ec3546b89998cfe44f45f1d0d9a875c73c/magtense-2.1.0+cuda116-py39-none-win_amd64.whl

SET wheels[3].name="2.1.0+cuda11.6_py3.10"
SET wheels[3].link=https://files.pythonhosted.org/packages/4c/a1/cab6e91b31e6642a86e0af24240e3d3d170ff65f549a78925bb6e0fb7e27/magtense-2.1.0+cuda116-py310-none-win_amd64.whl

SET wheels[4].name="2.1.0+cuda11.7_py3.9"
SET wheels[4].link=https://files.pythonhosted.org/packages/d6/77/5e19ad9683567d54bc02a6c50be2389b4135ae22d586bb66b0229c7574c2/magtense-2.1.0+cuda117-py39-none-win_amd64.whl

SET wheels[5].name="2.1.0+cuda11.7_py3.10"
SET wheels[5].link=https://files.pythonhosted.org/packages/f5/d8/c00e6973dd92ce6e903ae091cf03d040edf2e0767d93b92374d4ac89489f/magtense-2.1.0+cuda117-py310-none-win_amd64.whl

REM Installation of matching wheel
FOR /L %%i IN ( 0 1 5 ) DO  (
    IF !wheels[%%i].name! == "%MT_VERSION%-%CUDA_PIP%_py%PY_PIP%" (
        CALL python -m pip install -vv --use-deprecated=legacy-resolver %%wheels[%%i].link%%
    )
)
