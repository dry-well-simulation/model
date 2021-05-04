import numpy as np


def compute(dt, Qw, A, K_h, K_v, hw, R, r_w, zspan, L_w, n, theta_i, Z, PSI_i, r_initial, eps_R, eps_h, RK):
    """
    compute water level in well, position of wetting front (around and below the well), and spilled
    discharge rate at a single time step.

    :param dt: (float) the present time step
    :param Qw: (float) the discharge into the well
    :param A: (float) the well cross section area
    :param K_h: (float) the horizontal conductivity
    :param K_v: (float) the vertical conductivity
    :param hw: (float) the head at the well at the start of the time step
    :param R: (ndarray) the position of wetting front at the start of the time step
    :param r_w: (float) the radius of the well
    :param zspan: (ndarray) the elevation of the vertical numerical sections of the well
    :param L_w: (float) the elevation of the well
    :param n: (float) the porosity of the porous media
    :param theta_i: (float) the initial water content (theta_i<=n)
    :param Z: (float ??) the depth of the wetting front below the well, at the start of the time step, in the deepest point??
    :param PSI_i: the initial matric head
    :param r_initial: todo....
    :param eps_R: wetting front location tolerance [-]
    :param eps_h:  well water level tolerance [-]
    :param RK: number of Runge–Kutta sections [-]
    :return:
        next_hw : (float) hw at the end of the time step
        next_R: (ndarray) R at the end of the time step
        next_Q_spill: (float) The discharge rate of spill from the well top
        next_Z: (float) Z at the end of the time step
        next_zR (float) R at Z at the end of the time step
        dt2 (float) the time step used in practice
    """
    # initialize
    O_r = [1, 1]  # relative change of R between iterations... todo: find better name
    O_h = [1, 1]
    next_R = R
    next_hw = hw
    while O_r[0] > eps_R or O_h[0] > eps_h:
        temp_R = next_R
        temp_hw = next_hw
        O_r, next_R = next_radius(K_h, K_v, next_hw, hw, zspan, next_R, R, r_w, r_initial, PSI, n, theta_i, dt,
                                  temp_R, O_r, eps_R, RK)

        O_h, next_hw, next_Q_spill, next_Z, next_zR = \
            next_height(K_v, n, theta_i, dt, Z, next_R, R, zspan, r_w, hw, Qw, A, L_w, temp_hw, O_h)
        if O_h[1] >= O_h[0] or next_hw < 0:
            dt = dt / 2
            O_h[0] = 1
        else:
            O_h[0] = O_h[1]


    return None
    # return next_hw,next_R,next_Q_spill,next_Z,next_zR,dt2
