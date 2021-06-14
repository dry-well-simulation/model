from numpy import pi
import numpy as np


EPSILON = 0.01
EXAGGERATION_FACTOR = 10


def core_simulation(horiz_cond, well_length, simulation_time, well_radius, porosity, initial_VWC, recharge_duration,
                    recharge_rate, anisotropy_ratio=1, tol_radius=0.01, tol_head=0.01,
                    n_vertical_slices=100, RungeKutta_sections=4, initial_matric_head=0, output_times=None, dt=None):
    """
    Perform the entire calculation of the dry well recharge, for specified parameters.
    The units are written here as [m] and [min], but in practice they can be whatever
    length/ time units, as long as you are consistent.
    :param horiz_cond: float, horizontal hydraulic conductivity [m/min]
    :param well_length: float, the length of the well, assumed to be screened throughout  [m]
    :param simulation_time: float, total simulation time [min]
    :param well_radius: float, well's radius [m]
    :param porosity:  float, porosity [-]
    :param initial_VWC: float, initial volumetric water content [-]
    :param recharge_duration: float, time in which injection is on [min] # todo: make it possible to input nd array
    :param recharge_rate: discharge into the well when the discharge is on [m^3/min]
                        # todo: make it possible to input nd array
    :param anisotropy_ratio: float, the ratio of horizontal to vertical hydraulic conductivity [-]
    :param tol_radius: float, wetting front location tolerance [-]
    :param tol_head: float,  well water level tolerance [-]
    :param dt: float, length of time steps during simulation [-]
    :param n_vertical_slices: float, number of vertical slices along well [-]
    :param RungeKutta_sections: int, number of Runge–Kutta sections [-]
    :param initial_matric_head: float, initial matric head [m]. Must be negative.
    :param output_times:  float or nd array, the time(s) at which output is saved and returned [min]
                          if None, will use a default set.
    :return:
        R_plot, Z_plot: radius and vertical location of wetting front at times to plot [m]
        zspan: the vertical slices locations [m]
        hw_all: head at well, for any time in t [m]
        t: times at which hw_all is given [min]
    """
    assert horiz_cond > 0
    assert well_length > 0
    assert isinstance(n_vertical_slices, int)
    assert initial_matric_head <= 0
    vert_cond = horiz_cond / anisotropy_ratio  # vertical conductivity
    well_section_area = pi * well_radius ** 2
    r_initial = 1 + EPSILON  + (- initial_matric_head) / (well_length - initial_matric_head)
    # The initial relative wetting radius is defined such that r_initial > 1, but not much larger than 1.
    if dt is None:
        dt = well_section_area / \
             (EXAGGERATION_FACTOR * 2 * pi * horiz_cond * (well_length + - initial_matric_head))
    n_time_steps = int(np.ceil(simulation_time / dt))
    t = np.linspace(start=0, stop=simulation_time, num=n_time_steps)
    dt_all = np.diff(t)
    dz_max = recharge_duration * recharge_rate \
             / (EXAGGERATION_FACTOR * 2 * pi * well_radius * (well_length + - initial_matric_head))
    if well_length / n_vertical_slices > dz_max:
        n_vertical_slices = int(np.ceil(well_length / dz_max))
    zspan = np.linspace(0, well_length, n_vertical_slices)  # TODO: chg name

    Qw_all = np.where(t <= recharge_duration, recharge_rate, 0)
    hw_all = np.zeros((n_time_steps,))  # TODO: change to np.nan ?
    Q_spill = np.zeros((n_time_steps,))  # TODO: change to np.nan ?
    V_total = np.zeros((n_time_steps,))  # TODO: change to np.nan ?
    V_pm = np.zeros((n_time_steps,))  # TODO: change to np.nan ?
    R_all = np.zeros((n_vertical_slices, n_time_steps))  # TODO: change to np.nan ?
    R_all[:, 0] = well_radius
    Z_vertical = np.zeros((n_time_steps,))  # TODO: change to np.nan ?
    R_vertical = well_radius * np.ones((n_time_steps,))  # TODO: change to np.nan ?

    for ti in range(n_time_steps - 1):
        dt = dt_all[ti]
        Qw = (Qw_all[ti] + Qw_all[ti+1]) / 2  # present Qw
        hw = hw_all[ti]
        R = R_all[:, ti]
        Z = Z_vertical[ti]

        next_hw, next_R, next_Q_spill, next_Z, next_zR, dt2 = \
            compute(dt, Qw, well_section_area, horiz_cond, vert_cond, hw, R, well_radius, zspan, well_length, porosity, initial_VWC, Z, initial_matric_head,
                    r_initial, tol_radius, tol_head, RungeKutta_sections)
        # todo: chg to "compute_timestep"
        if dt2 < dt:
            dt2_2 = dt2
            while dt2_2 <= dt:
                next_hw, next_R, next_Q_spill, next_Z, next_zR, dt2 = \
                    compute(dt2, Qw, well_section_area, horiz_cond, vert_cond, next_hw, next_R, well_radius, zspan, well_length, porosity, initial_VWC, next_Z, initial_matric_head,
                            r_initial, tol_radius, tol_head, RungeKutta_sections)
                dt2_2 = dt2_2 + dt2
        hw_all[ti + 1] = next_hw
        R_all[:, ti + 1] = next_R
        Q_spill[ti + 1] = next_Q_spill
        Z_vertical[ti + 1] = next_Z
        R_vertical[ti + 1] = next_zR
        V_total[ti + 1] = V_total[ti] + Qw * dt  # total volume of water that was recharged into the well
        V_pm[ti + 1] = V_total[ti + 1] - next_hw * well_section_area # water volume in porous media is total vol. minus vol. in well.

    Z_vert_origin = Z_vertical - max(Z_vertical)  # construct the vertical flow as a reconstruction of the bottom R ????
    # R_v = R_all[0, :]

    if output_times is None:
        output_times = [simulation_time]  # times to show in figure
        if recharge_duration < simulation_time:
            output_times.append(recharge_duration)
        output_times = np.sort(output_times)

    # h_w = hw_all
    Z_vert_fl = np.flip(Z_vert_origin)
    R_vert_fl = well_radius * np.ones((n_time_steps, len(output_times)))
    R_plot = dict()
    Z_plot = dict()
    for i, t_to_show_now in enumerate(output_times):
        idx_of_t_to_show_now = np.argwhere(t >= t_to_show_now)[0][0]  # todo: what if this fails? do something
        R_vert_fl_tmp = np.flip(R_vertical[:idx_of_t_to_show_now])
        R_vert_fl[:idx_of_t_to_show_now, i] = R_vert_fl_tmp  # not sure this will work
        R_plot[t_to_show_now] = np.concatenate((R_vert_fl[:, i], R_all[:, idx_of_t_to_show_now],))
        Z_plot[t_to_show_now] = np.concatenate((Z_vert_fl.T, zspan, ))
    return R_plot, Z_plot, zspan, hw_all, t


def test_core_simulation():
    R_plot, Z_plot, zspan, hw_all, t = \
        core_simulation(horiz_cond=5 / 24 / 60, well_length=10, simulation_time=400, well_radius=0.2, porosity=0.3,
                        initial_VWC=0.1, recharge_duration=30, recharge_rate=30 / 60)




def calculate_Q_slice(K, hw, zspan, R, r_w, r_initial, PSI):
    R_correct = R
    Q_slice = np.zeros_like(zspan)
    z_below_water = zspan < hw
    if np.any(z_below_water):
        z_below_water[np.argwhere(z_below_water)[-1][0]] = False

    for i, z in enumerate(zspan[z_below_water]):
        eps_x = R[i] / r_w
        if eps_x > r_initial:
            Q_slice[i] = 2 * pi * K * (hw - z - PSI) / np.log(eps_x)
        else:
            eps_x = r_initial
            R_correct[i] = r_w * eps_x
            Q_slice[i] = 2 * pi * K * (hw - z - PSI) / np.log(eps_x)
        if Q_slice[i] < 0:
            print('never get her. todo remov')
            Q_slice[i] = 0
    return Q_slice, R_correct


def calculate_Q_v_slice(K_v, R, r_w):
    # calc the net vertical flow from each slice, acounting for the slice above
    Q_v_slice_all = np.zeros_like(R)
    R_tmp = R.T  # todo: do we need this?
    calc_at = R_tmp > r_w
    Q_v_slice_all[calc_at] = -pi * (R_tmp[calc_at] ** 2 - r_w ** 2) * K_v
    single_zero = np.array([0])
    return np.concatenate((-np.diff(Q_v_slice_all), single_zero))


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


def next_height(K_v, n, theta_i, dt, Z, next_R, R, zspan, r_w, hw, Qw, A, L_w, temp_hw, O_h):
    dz = K_v / (n - theta_i) * dt
    next_Z = Z + dz
    next_zR = next_R[0]
    Q_pm = np.sum(pi * (next_R ** 2 - R ** 2) * (n - theta_i) * zspan[1]) + dz * pi * (next_zR ** 2 - r_w ** 2) * (n - theta_i)
    if Q_pm < 0:
        Q_pm = 0
    next_hw_theoretical = hw + (Qw * dt - Q_pm) / A
    if next_hw_theoretical > L_w:
        next_Q_spill = (next_hw_theoretical - L_w) * A / dt
        next_hw = L_w
    else:
        next_Q_spill = 0
        next_hw = next_hw_theoretical
    O_h[1] = np.abs(next_hw - temp_hw) / L_w

    return O_h, next_hw,  next_Q_spill, next_Z, next_zR


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
