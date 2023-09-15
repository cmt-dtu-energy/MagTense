#!/bin/bash

# Define dictornary of wheels on PyPI
declare -A wheels

wheels['2.1.0+cpu_py3.9_unix']='https://files.pythonhosted.org/packages/1f/9d/682557bbaea87ceaf2b21a6cdc345e7e2b108d711f07c5228a94c6b31295/magtense-2.1.0+cpu-py39-none-manylinux1_x86_64.whl'
wheels['2.1.0+cpu_py3.9_win']='https://files.pythonhosted.org/packages/56/1a/e7bde7359e07318d610d05ba52484fb5e7318ce693360c1ad6ab26342e22/magtense-2.1.0+cpu-py39-none-win_amd64.whl'
wheels['2.1.0+cpu_py3.10_unix']='https://files.pythonhosted.org/packages/7f/df/278960701ef5a0c78cd4a44eb264c230a8d27206265b470cd8af89409bad/magtense-2.1.0+cpu-py310-none-manylinux1_x86_64.whl'
wheels['2.1.0+cpu_py3.10_win']='https://files.pythonhosted.org/packages/7a/f6/df4567440d519873c780365848e9beca2f59e71ae0818953d872a928c2b5/magtense-2.1.0+cpu-py310-none-win_amd64.whl'

wheels['2.1.0+cuda11.4_py3.9_unix']='https://files.pythonhosted.org/packages/06/c2/5dff71989e6ad559073535067888d0096a41fc497b699fe4ec9b7abb4428/magtense-2.1.0+cuda114-py39-none-manylinux1_x86_64.whl'
wheels['2.1.0+cuda11.4_py3.10_unix']='https://files.pythonhosted.org/packages/25/7b/39a45d45953d86a50c86ba0d1690b4ea9c3836196535f9615db62ac850e4/magtense-2.1.0+cuda114-py310-none-manylinux1_x86_64.whl'

wheels['2.1.0+cuda11.6_py3.9_unix']='https://files.pythonhosted.org/packages/d3/c5/c9f75b68d108a7d07e3fd3213681b3dad778f5980f3409319d87661fcf4d/magtense-2.1.0+cuda116-py39-none-manylinux1_x86_64.whl'
wheels['2.1.0+cuda11.6_py3.10_unix']='https://files.pythonhosted.org/packages/c5/aa/7c6056a0bda349fbc90ecea115fad7822b8050839ee3a5938f6ecc55f7bb/magtense-2.1.0+cuda116-py310-none-manylinux1_x86_64.whl'
wheels['2.1.0+cuda11.6_py3.9_win']='https://files.pythonhosted.org/packages/ab/06/7c1b285066e21a0ef19fedb540ec3546b89998cfe44f45f1d0d9a875c73c/magtense-2.1.0+cuda116-py39-none-win_amd64.whl'
wheels['2.1.0+cuda11.6_py3.10_win']='https://files.pythonhosted.org/packages/4c/a1/cab6e91b31e6642a86e0af24240e3d3d170ff65f549a78925bb6e0fb7e27/magtense-2.1.0+cuda116-py310-none-win_amd64.whl'

wheels['2.1.0+cuda11.7_py3.9_unix']='https://files.pythonhosted.org/packages/6e/ec/caa8adcb3ff7829ed9d1cfed765901e4c4d5cabf57c3ce37b65a0755d7f4/magtense-2.1.0+cuda117-py39-none-manylinux1_x86_64.whl'
wheels['2.1.0+cuda11.7_py3.10_unix']='https://files.pythonhosted.org/packages/2f/bc/4deb21e25c9a7b319b432ae20873852fb73293f16f5ea2efb2532b3ccc9d/magtense-2.1.0+cuda117-py310-none-manylinux1_x86_64.whl'
wheels['2.1.0+cuda11.7_py3.9_win']='https://files.pythonhosted.org/packages/d6/77/5e19ad9683567d54bc02a6c50be2389b4135ae22d586bb66b0229c7574c2/magtense-2.1.0+cuda117-py39-none-win_amd64.whl'
wheels['2.1.0+cuda11.7_py3.10_win']='https://files.pythonhosted.org/packages/f5/d8/c00e6973dd92ce6e903ae091cf03d040edf2e0767d93b92374d4ac89489f/magtense-2.1.0+cuda117-py310-none-win_amd64.whl'

wheels['2.2.0-cpu_py3.9_unix']='https://files.pythonhosted.org/packages/27/c7/e14a412c0400e755a70b0cacd324f6050918b57b99f21a9e5850908c2477/magtense-2.2.0-0-py39-none-manylinux1_x86_64.whl'
wheels['2.2.0-cpu_py3.10_unix']='https://files.pythonhosted.org/packages/8d/d5/db1845cf9d43f6cb8ea9338ee13fe2de553370a388bdec65fb8be1ea0126/magtense-2.2.0-0-py310-none-manylinux1_x86_64.whl'
wheels['2.2.0-cpu_py3.11_unix']='https://files.pythonhosted.org/packages/da/89/6a2ab09bfd31576d72d737b44c7ab098cd44a51b237a06d501e2b3a52585/magtense-2.2.0-0-py311-none-manylinux1_x86_64.whl'
wheels['2.2.0-cpu_py3.9_win']='https://files.pythonhosted.org/packages/6e/f2/b8203604efe834d6dc1b2634c4fb177557fdc27a0556f3317bdeced3d63a/magtense-2.2.0-0-py39-none-win_amd64.whl'
wheels['2.2.0-cpu_py3.10_win']='https://files.pythonhosted.org/packages/c9/a6/763b9eea716f9969b65fdfe938d87b89245f5b42dc0c2a6462695f40f03e/magtense-2.2.0-0-py310-none-win_amd64.whl'
wheels['2.2.0-cpu_py3.11_win']='https://files.pythonhosted.org/packages/d2/14/55cbf1fc9222a6dd0df0a84317cb7af818bf38d4c5bf83d688b73417bc07/magtense-2.2.0-0-py311-none-win_amd64.whl'

wheels['2.2.0-cuda11_py3.9_unix']='https://files.pythonhosted.org/packages/6b/d4/95c58c6520c31b37af675ec69e1aa6234f51b9f2492fe86cff5a79470eff/magtense-2.2.0-1-py39-none-manylinux1_x86_64.whl'
wheels['2.2.0-cuda11_py3.10_unix']='https://files.pythonhosted.org/packages/37/ed/bcd9b48ae7aab9b87f0f1409c228aa1c7ed0a8864d75323f8dc60ad4dfbb/magtense-2.2.0-1-py310-none-manylinux1_x86_64.whl'
wheels['2.2.0-cuda11_py3.11_unix']='https://files.pythonhosted.org/packages/18/c1/8c2264c30e016e5df75515c3e0238d0cf9c94c056e268335d69f44d68853/magtense-2.2.0-1-py311-none-manylinux1_x86_64.whl'
wheels['2.2.0-cuda11_py3.9_win']='https://files.pythonhosted.org/packages/a0/56/a2dae121b517a4bdb741c0bfe74957c53c703e91f674603c1d8ec1c4c1cf/magtense-2.2.0-1-py39-none-win_amd64.whl'
wheels['2.2.0-cuda11_py3.10_win']='https://files.pythonhosted.org/packages/e9/f1/2c370a3f5d7ef381a2cd3adadc24256ab6d0e1ba4c493658a5fe5f2a3c8d/magtense-2.2.0-1-py310-none-win_amd64.whl'
wheels['2.2.0-cuda11_py3.11_win']='https://files.pythonhosted.org/packages/50/b0/1fdc988e23af9443de1a9dd40509d8bd04fd8e527dce96fbc7c5fd94a881/magtense-2.2.0-1-py311-none-win_amd64.whl'

wheels['2.2.0-cuda12_py3.9_unix']='https://files.pythonhosted.org/packages/b5/84/2956391ff6a08ae702d7ccacebc332dd32f266153340142b02aa1d922199/magtense-2.2.0-2-py39-none-manylinux1_x86_64.whl'
wheels['2.2.0-cuda12_py3.10_unix']='https://files.pythonhosted.org/packages/55/5e/f4267aeb5ad339c850c9338a04f5f586c8bbf769fcbd8ee3c66a477321a1/magtense-2.2.0-2-py310-none-manylinux1_x86_64.whl'
wheels['2.2.0-cuda12_py3.11_unix']='https://files.pythonhosted.org/packages/52/bb/d5bba58cf6393378261bbe07f8c2d9a80e18493cc0a99b3e4eed844dadd4/magtense-2.2.0-2-py311-none-manylinux1_x86_64.whl'
wheels['2.2.0-cuda12_py3.9_win']='https://files.pythonhosted.org/packages/79/27/f9b3c391ca48f82af03389bab68a691ad98fae7b81e3f424d38a3d2b7910/magtense-2.2.0-2-py39-none-win_amd64.whl'
wheels['2.2.0-cuda12_py3.10_win']='https://files.pythonhosted.org/packages/55/5e/f4267aeb5ad339c850c9338a04f5f586c8bbf769fcbd8ee3c66a477321a1/magtense-2.2.0-2-py310-none-manylinux1_x86_64.whl'
wheels['2.2.0-cuda12_py3.11_win']='https://files.pythonhosted.org/packages/92/b4/7752c40a40edfc84dcd673855b91c8e64f4b3e778e9a0604d313238edf3d/magtense-2.2.0-2-py311-none-win_amd64.whl'

# Installation of matching wheel
python -m pip install -vv --use-deprecated=legacy-resolver "${wheels[${MT_VERSION}-${CUDA_PIP}_py${PY_PIP}_${ARCH_PIP}]}"
