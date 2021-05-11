from numpy import pi
import numpy as np


def calculate_Q_slice(K_h, hw, zspan, R, r_w, r_initial, PSI):
    R_correct = R
    Q_slice = np.zeros_like(zspan)
    z_below_water = zspan < hw

    if np.any(z_below_water):
        # stop here
        if length(idx) > 1
            idx = idx(1:end - 1);
            else
            idx = [];

    return Q_slice, R_correct


def next_radius(K_h, K_v, next_hw, hw, zspan, next_R, R, r_w, r_initial, PSI, n, theta_i, dt, temp_R, O_r, eps_R, RK):
    # todo: use dz instead of zspan[1]
    Q_h_slice, next_R = calculate_Q_slice(K_h, (next_hw + hw) / 2, zspan, (next_R + R)/ 2, r_w, r_initial, PSI)
    Q_v_slice = calculate_Q_v_slice(K_v, (next_R + R) / 2, r_w).T
    next_R = R + (Q_h_slice + Q_v_slice / zspan[1]).T / (2 * pi * ((next_R + R) / 2) * (n-theta_i)) * dt
    O_r[1] = np.mean(np.abs(next_R - temp_R) / R)
    if O_r[1] < O_r[0]:
        O_r[0] = O_r[1]
    else:
        while O_r[1] > eps_R and O_r[1] > O_r[0]:
            temp_R = next_R
            R_T = R
            next_R_Q = R + (Q_h_slice + Q_v_slice/ zspan[1]).T / (2 * pi * R * (n-theta_i)) * (dt / RK)
            for _ in range(RK - 2):
                R_T = next_R_Q
                next_R_Q = R_T + (Q_h_slice + Q_v_slice / zspan[1]).T / (2 * pi * ((R_T + next_R_Q) / 2) * (n-theta_i)) * (dt / RK)

            next_R = next_R_Q + (Q_h_slice + Q_v_slice / zspan[1]).T / (2 * pi * ((next_R + R_T) / 2) * (n-theta_i)) * (dt / RK)
            O_r[1] = np.mean(np.abs(next_R - temp_R) / R)
        Q_h_slice, next_R = calculate_Q_slice(K_h, (next_hw + hw) / 2, zspan, (next_R + R) / 2, r_w, r_initial, PSI)
        Q_v_slice = calculate_Q_v_slice(K_v, (next_R + R) / 2, r_w).T
        next_R = R + (Q_h_slice - Q_v_slice / zspan[1]).T/ (2 * pi * ((next_R + R) / 2) * (n - theta_i)) * dt
        O_r[0] = np.mean(np.abs(next_R - temp_R) / R)
    return O_r, next_R

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
        O_r, next_R = next_radius(K_h, K_v, next_hw, hw, zspan, next_R, R, r_w, r_initial, PSI_i, n, theta_i, dt,
                                  temp_R, O_r, eps_R, RK)

        O_h, next_hw, next_Q_spill, next_Z, next_zR = \
            next_height(K_v, n, theta_i, dt, Z, next_R, R, zspan, r_w, hw, Qw, A, L_w, temp_hw, O_h)
        if O_h[1] >= O_h[0] or next_hw < 0:
            dt = dt / 2
            O_h[0] = 1
        else:
            O_h[0] = O_h[1]


    return next_hw, next_R, next_Q_spill, next_Z, next_zR, dt
