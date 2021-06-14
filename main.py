from scripts import core_simulation
import matplotlib.pyplot as plt


if __name__ == '__main__':
    # example run
    R_plot, Z_plot, zspan, hw_all, t = \
        core_simulation(horiz_cond=5 / 24 / 60, well_length=10, simulation_time=400, well_radius=0.2, porosity=0.3,
                        initial_VWC=0.1, recharge_duration=30, recharge_rate=30 / 60)

    plt.plot(t, hw_all, '.', label='h_w(t)')
    plt.legend()
    plt.show()

    # for time_ in R_plot.keys():
    #     plt.scatter(R_plot[time_], Z_plot[time_], label=time_)
    # plt.legend()
    # plt.show()

