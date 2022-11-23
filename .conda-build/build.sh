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

# Installation of matching wheel
python -m pip install -vv --use-deprecated=legacy-resolver "${wheels[${MT_VERSION}+${CUDA_PIP}_py${PY_PIP}_${ARCH_PIP}]}"
