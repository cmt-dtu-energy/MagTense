def compute_magnetic_energy(M, H_exc, H_ext, H_dem, H_ani, Vols, problem, Ms):
    import numpy as np

    mu0 = 4 * np.pi * 1e-7
    m = np.zeros_like(M)

    # loop over first, third and fourth dimensions of m
    for i in range(m.shape[0]):
        for k in range(m.shape[2]):
            for l in range(m.shape[3]):
                m[i, :, k, l] = mu0 * Ms * Vols * M[i, :, k, l]

    # elementwise dot-product m . H summed over components, then sum over spatial cells
    E_exc = 0.5 * np.sum(m[:, :, :, 0] * H_exc[:, :, :, 0] + m[:, :, :, 1] * H_exc[:, :, :, 1] + m[:, :, :, 2] * H_exc[:, :, :, 2], axis=1)  # [J]
    E_ext =       np.sum(m[:, :, :, 0] * H_ext[:, :, :, 0] + m[:, :, :, 1] * H_ext[:, :, :, 1] + m[:, :, :, 2] * H_ext[:, :, :, 2], axis=1)       # [J]
    E_dem = 0.5 * np.sum(m[:, :, :, 0] * H_dem[:, :, :, 0] + m[:, :, :, 1] * H_dem[:, :, :, 1] + m[:, :, :, 2] * H_dem[:, :, :, 2], axis=1) # [J]
    E_ani = 0.5 * np.sum(m[:, :, :, 0] * H_ani[:, :, :, 0] + m[:, :, :, 1] * H_ani[:, :, :, 1] + m[:, :, :, 2] * H_ani[:, :, :, 2], axis=1) # [J]

    # add constant anisotropy K0 * Vols
    E_ani = E_ani + np.sum(((np.asarray(problem.K0)).ravel() * np.asarray(Vols).ravel()).reshape(-1, 1))

    # energy densities (divide total energy by total volume)
    total_volume = np.prod(np.asarray(problem.grid_L))
    E_exc_dens = E_exc / total_volume
    E_ext_dens = E_ext / total_volume
    E_dem_dens = E_dem / total_volume
    E_ani_dens = E_ani / total_volume

    # characteristic energy density Km
    Km = 0.5 * (Ms ** 2) * mu0

    E_exc_red = E_exc_dens / Km
    E_ext_red = E_ext_dens / Km
    E_dem_red = E_dem_dens / Km
    E_ani_red = E_ani_dens / Km

    return E_dem_red, E_exc_red, E_ani_red, E_ext_red