import numpy as np
import pysixtrack


def vectorize_all_coords(Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                         Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                         Dsigma_wrt_CO_m, Ddelta_wrt_CO):

    n_part = 1
    for var in [Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                Dsigma_wrt_CO_m, Ddelta_wrt_CO]:
        if hasattr(var, '__iter__'):
            if n_part == 1:
                n_part = len(var)
            assert len(var) == n_part

    Dx_wrt_CO_m = Dx_wrt_CO_m + np.zeros(n_part)
    Dpx_wrt_CO_rad = Dpx_wrt_CO_rad + np.zeros(n_part)
    Dy_wrt_CO_m = Dy_wrt_CO_m + np.zeros(n_part)
    Dpy_wrt_CO_rad = Dpy_wrt_CO_rad + np.zeros(n_part)
    Dsigma_wrt_CO_m = Dsigma_wrt_CO_m + np.zeros(n_part)
    Ddelta_wrt_CO = Ddelta_wrt_CO + np.zeros(n_part)

    return Dx_wrt_CO_m, Dpx_wrt_CO_rad,\
     Dy_wrt_CO_m, Dpy_wrt_CO_rad,\
     Dsigma_wrt_CO_m, Ddelta_wrt_CO



def track_particle_pysixtrack(line, part, Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                              Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                              Dsigma_wrt_CO_m, Ddelta_wrt_CO, n_turns, verbose=False):

    Dx_wrt_CO_m, Dpx_wrt_CO_rad,\
        Dy_wrt_CO_m, Dpy_wrt_CO_rad,\
        Dsigma_wrt_CO_m, Ddelta_wrt_CO = vectorize_all_coords(
                             Dx_wrt_CO_m, Dpx_wrt_CO_rad,
                             Dy_wrt_CO_m, Dpy_wrt_CO_rad,
                             Dsigma_wrt_CO_m, Ddelta_wrt_CO)

    part.x += Dx_wrt_CO_m
    part.px += Dpx_wrt_CO_rad
    part.y += Dy_wrt_CO_m
    part.py += Dpy_wrt_CO_rad
    part.sigma += Dsigma_wrt_CO_m
    part.delta += Ddelta_wrt_CO

    x_tbt = []
    px_tbt = []
    y_tbt = []
    py_tbt = []
    sigma_tbt = []
    delta_tbt = []

    for i_turn in range(n_turns):
        if verbose:
            print('Turn %d/%d' % (i_turn, n_turns))

        x_tbt.append(part.x.copy())
        px_tbt.append(part.px.copy())
        y_tbt.append(part.y.copy())
        py_tbt.append(part.py.copy())
        sigma_tbt.append(part.sigma.copy())
        delta_tbt.append(part.delta.copy())

        line.track(part)

    x_tbt = np.array(x_tbt)
    px_tbt = np.array(px_tbt)
    y_tbt = np.array(y_tbt)
    py_tbt = np.array(py_tbt)
    sigma_tbt = np.array(sigma_tbt)
    delta_tbt = np.array(delta_tbt)

    return x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt
