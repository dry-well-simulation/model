from scripts import core_simulation
import matplotlib.pyplot as plt


if __name__ == '__main__':
    # example run
    R_plot, Z_plot, zspan, hw_all, t = \
        core_simulation(K=5/24/60, L_w=10, T=400, r_w=0.2, n=0.3,
                        theta_i=0.1, T_on=30, Qw_on=30/60, Nt=1000)

    plt.plot(t, hw_all, '.', label='h_w(t)')
    plt.legend()
    plt.show()

    for time_ in R_plot.keys():
        plt.scatter(R_plot[time_], Z_plot[time_], label=time_)
    plt.legend()
    plt.show()

