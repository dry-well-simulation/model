import numpy as np
from scripts import compute
import matplotlib.pyplot as plt


EPSILON = 0.01
EXAGGERATION_FACTOR = 10


def main(K, L_w, T, r_w, n, theta_i, T_on, Qw_on, I=1, eps_R=0.01, eps_h=0.01, time_int=1000, Nz=100, RK=4, PSI_i=0, t_to_show=None):
    """ example ziv
    perform the entire calculation of the well recharge
    :param K: horizontal hydraulic conductivity [m/min]
    :param L_w: well's length [m]
    :param T: total simulation time [min]
    :param r_w: well's radius [m]
    :param n:  porosity [-]
    :param theta_i: initial water content [-]
    :param T_on: time in which injection is on [min]
    :param Qw_on: discharge into the well when the discharge is on [m^3/min]
    :param I:  anisotropy ratio [-]
    :param eps_R: wetting front location tolerance [-]
    :param eps_h:  well water level tolerance [-]
    :param time_int: number of time intervals [-]
    :param Nz: number of vertical slices along well [-]
    :param RK: number of Runge–Kutta sections [-]
    :param PSI_i: initial matric head [m]
    :return:
        R_plot, Z_plot: radius and vertical location of wetting front at times to plot [m]
        zspan: the vertical slices locations [m]
        hw_all: head at well, for any time in t [m]
        t: times at which hw_all is given [min]
    """
    assert K > 0
    assert L_w > 0
    assert isinstance(Nz, int)
    assert PSI_i <= 0
    K_h = K  # horizontal conductivity
    K_v = K / I  # vertical conductivity
    A = np.pi * r_w ** 2
    r_initial = (1 + EPSILON) - PSI_i / (L_w - PSI_i)
    # initial relative wetting radius (r_initial > 1, but not much larger than 1)
    dt_max = A / (EXAGGERATION_FACTOR * 2 * np.pi * K_h * (L_w - PSI_i))
    if T/time_int > dt_max:
        time_int = int(np.ceil(T / dt_max))  #TODO: call time_int with different name
    t = np.linspace(0, T, time_int)
    dt_all = np.diff(t)
    dz_max = T_on * Qw_on / (EXAGGERATION_FACTOR * 2 * np.pi * r_w * (L_w - PSI_i))
    if L_w / Nz > dz_max:
        Nz = int(np.ceil(L_w / dz_max))
    zspan = np.linspace(0, L_w, Nz)  # TODO: chg name

    Qw_all = np.where(t <= T_on, Qw_on, 0)
    hw_all = np.zeros((time_int,))  # TODO: change to np.nan ?
    Q_spill = np.zeros((time_int,))  # TODO: change to np.nan ?
    V_total = np.zeros((time_int,))  # TODO: change to np.nan ?
    V_pm = np.zeros((time_int,))  # TODO: change to np.nan ?
    R_all = np.zeros((Nz, time_int))  # TODO: change to np.nan ?
    R_all[:, 0] = r_w
    Z_vertical = np.zeros((time_int,))  # TODO: change to np.nan ?
    R_vertical = r_w * np.ones((time_int,))  # TODO: change to np.nan ?

    for ti in range(time_int-1):
        dt = dt_all[ti]
        Qw = (Qw_all[ti] + Qw_all[ti+1]) / 2  # present Qw
        hw = hw_all[ti]
        R = R_all[:, ti]
        Z = Z_vertical[ti]

        next_hw, next_R, next_Q_spill, next_Z, next_zR, dt2 = \
            compute(dt, Qw, A, K_h, K_v, hw, R, r_w, zspan, L_w, n, theta_i, Z, PSI_i,
                    r_initial, eps_R, eps_h, RK)
        # todo: chg to "compute_timestep"
        if dt2 < dt:
            dt2_2 = dt2
            while dt2_2 <= dt:
                next_hw, next_R, next_Q_spill, next_Z, next_zR, dt2 = \
                    compute(dt2, Qw, A, K_h, K_v, next_hw, next_R, r_w, zspan, L_w, n, theta_i, next_Z, PSI_i,
                            r_initial, eps_R, eps_h,RK)
                dt2_2 = dt2_2 + dt2
        hw_all[ti + 1] = next_hw
        R_all[:, ti + 1] = next_R
        Q_spill[ti + 1] = next_Q_spill
        Z_vertical[ti + 1] = next_Z
        R_vertical[ti + 1] = next_zR
        V_total[ti + 1] = V_total[ti] + Qw * dt # total volume of water that was recharged into the well
        V_pm[ti + 1] = V_total[ti + 1] - next_hw * A # water volume in porous media is total vol. minus vol. in well.

    Z_vert_origin = Z_vertical - max(Z_vertical) # construct the vertical flow as a reconstruction of the bottom R ????
    # R_v = R_all[0, :]

    if t_to_show is None:
        t_to_show = np.array([33, 90, 120, 180, 360]) #  times to show in figure
        # todo: make it relative to T_on, not abs values
    # h_w = hw_all
    Z_vert_fl = np.flip(Z_vert_origin)
    R_vert_fl = r_w * np.ones((time_int, len(t_to_show)))
    R_plot = dict()
    Z_plot = dict()
    for i, t_to_show_now in enumerate(t_to_show):
        idx_of_t_to_show_now = np.argwhere(t >= t_to_show_now)[0][0]  # todo: what if this fails? do something
        R_vert_fl_tmp = np.flip(R_vertical[:idx_of_t_to_show_now])
        R_vert_fl[:idx_of_t_to_show_now, i] = R_vert_fl_tmp  # not sure this will work
        R_plot[t_to_show_now] = np.concatenate((R_vert_fl[:, i], R_all[:, idx_of_t_to_show_now],))
        Z_plot[t_to_show_now] = np.concatenate((Z_vert_fl.T, zspan, ))
    return R_plot, Z_plot, zspan, hw_all, t


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    R_plot, Z_plot, zspan, hw_all, t = \
        main(K=5/24/60, L_w=10, T=400, r_w=0.2, n=0.3, theta_i=0.1, T_on=30, Qw_on=30/60)

    plt.plot(t, hw_all, '.', label='h_w(t)')
    plt.legend()
    plt.show()

    for time_ in R_plot.keys():
        plt.scatter(R_plot[time_], Z_plot[time_], label=time_)
    plt.legend()
    plt.show()

