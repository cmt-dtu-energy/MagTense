import numpy as np
import matplotlib.pyplot as plt
from magtense.magstatics import Tiles, run_simulation

def read_three_columns(filename):
    col1 = np.array([])
    col2 = np.array([])
    col3 = np.array([])
    try:
        with open(filename, 'r') as file:
            for line_number, line in enumerate(file, start=1):
                # Skip empty lines
                if line.strip() == "":
                    continue

                parts = line.strip().split()

                if len(parts) != 3:
                    #print(f"Skipping line {line_number}: expected 3 columns, got {len(parts)}")
                    continue

                try:
                    value1 = float(parts[0])
                    value2 = float(parts[1])
                    value3 = float(parts[2])
                    col1=np.append(col1,value1)
                    col2=np.append(col2,value2)
                    col3=np.append(col3,value3)
                except ValueError:
                    print(f"Skipping line {line_number}: non-numeric data '{parts}'")

        return col1, col2, col3
    except FileNotFoundError:
        print(f"File '{filename}' not found.")
        return [], []

mu0 = 4 * np.pi * 1e-7

# Observation points
x = np.linspace(-20, 20, 201)

pts = np.column_stack([
    x,
    np.full_like(x, 5.0),
    np.full_like(x, 5.0)
]).astype(np.float64, order='F')

obs_size = np.column_stack([
    np.full_like(x, 2.0),
    np.full_like(x, 4.0),
    np.full_like(x, 3.0)
]).astype(np.float64, order='F')

# Magnetization directions to test
easy_axes = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
]

data = []


for easy_axis in easy_axes:

    tiles = Tiles(
        n=1,
        M_rem=1,
        easy_axis=easy_axis,
        tile_type=[8],
        color=[[1, 0, 0]]
    )

    tiles.size = ([3.0, 5.0, 1.0], 0)
    tiles.offset = ([0.0, 0.0, 0.0], 0)

    updated_tiles, H_demag = run_simulation(
        tiles,
        pts,
        obs_size=obs_size
    )

    data.append(H_demag*mu0)

val_path = "../../../documentation/examples_FEM_validation/Validation_avgprism/"
data=np.array(data)
xs_comsol = read_three_columns(val_path+"Avg_validation_comsol.txt")[0]
xx_comsol = read_three_columns(val_path+"Avg_validation_comsol.txt")[1]
xy_comsol = read_three_columns(val_path+"Avg_validation_comsol.txt")[2]

diff_xx=abs(xx_comsol-data[0][:, 0])
diff_xy = abs(xy_comsol-data[0][:, 1])

print("xx_err= ",np.sum(diff_xx)," , xy_err = ", np.sum(diff_xy))
#print(data)
plt.figure(figsize=(8, 5))
plt.plot(xs_comsol,xx_comsol,'ko',label=r"$\mathcal{N}_{xx}$, COMSOL")
plt.plot(xs_comsol,xy_comsol,'ro',label=r"$\mathcal{N}_{xy}$, COMSOL")
plt.plot(x, data[0][:, 0],'k', label=r"$\langle \mathcal{N}\rangle_{xx}$")
plt.plot(x, data[0][:, 1],'r', label=r"$\langle \mathcal{N}\rangle_{xy}$")
#plt.plot(x, data[0][:, 2], label='XZ')
#plt.plot(x, data[1][:, 1], label='YY')
#plt.plot(x, data[2][:, 1], label='ZY')
#plt.plot(x, data[2][:, 2], label='ZZ')

plt.xlabel("x")
plt.ylabel("H_demag*mu0")

plt.legend()
plt.grid(True)
plt.show()