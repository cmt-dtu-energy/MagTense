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

SET wheels[6].name="2.2.0-cpu_py3.9"
SET wheels[6].link=https://files.pythonhosted.org/packages/6e/f2/b8203604efe834d6dc1b2634c4fb177557fdc27a0556f3317bdeced3d63a/magtense-2.2.0-0-py39-none-win_amd64.whl

SET wheels[7].name="2.2.0-cpu_py3.10"
SET wheels[7].link=https://files.pythonhosted.org/packages/c9/a6/763b9eea716f9969b65fdfe938d87b89245f5b42dc0c2a6462695f40f03e/magtense-2.2.0-0-py310-none-win_amd64.whl

SET wheels[8].name="2.2.0-cpu_py3.11"
SET wheels[8].link=https://files.pythonhosted.org/packages/d2/14/55cbf1fc9222a6dd0df0a84317cb7af818bf38d4c5bf83d688b73417bc07/magtense-2.2.0-0-py311-none-win_amd64.whl

SET wheels[9].name="2.2.0-cuda11_py3.9"
SET wheels[9].link=https://files.pythonhosted.org/packages/a0/56/a2dae121b517a4bdb741c0bfe74957c53c703e91f674603c1d8ec1c4c1cf/magtense-2.2.0-1-py39-none-win_amd64.whl

SET wheels[10].name="2.2.0-cuda11_py3.10"
SET wheels[10].link=https://files.pythonhosted.org/packages/e9/f1/2c370a3f5d7ef381a2cd3adadc24256ab6d0e1ba4c493658a5fe5f2a3c8d/magtense-2.2.0-1-py310-none-win_amd64.whl

SET wheels[11].name="2.2.0-cuda11_py3.11"
SET wheels[11].link=https://files.pythonhosted.org/packages/50/b0/1fdc988e23af9443de1a9dd40509d8bd04fd8e527dce96fbc7c5fd94a881/magtense-2.2.0-1-py311-none-win_amd64.whl

SET wheels[12].name="2.2.0-cuda12_py3.9"
SET wheels[12].link=https://files.pythonhosted.org/packages/79/27/f9b3c391ca48f82af03389bab68a691ad98fae7b81e3f424d38a3d2b7910/magtense-2.2.0-2-py39-none-win_amd64.whl

SET wheels[13].name="2.2.0-cuda12_py3.10"
SET wheels[13].link=https://files.pythonhosted.org/packages/ca/64/f0a53f99b2ce1c175d8ea602aaaa28c3494b4cfd383df7c341fe72da1155/magtense-2.2.0-2-py310-none-win_amd64.whl

SET wheels[14].name="2.2.0-cuda12_py3.11"
SET wheels[14].link=https://files.pythonhosted.org/packages/92/b4/7752c40a40edfc84dcd673855b91c8e64f4b3e778e9a0604d313238edf3d/magtense-2.2.0-2-py311-none-win_amd64.whl

REM Installation of matching wheel
FOR /L %%i IN ( 0 1 5 ) DO  (
    IF !wheels[%%i].name! == "%MT_VERSION%-%CUDA_PIP%_py%PY_PIP%" (
        CALL python -m pip install -vv --use-deprecated=legacy-resolver %%wheels[%%i].link%%
    )
)
