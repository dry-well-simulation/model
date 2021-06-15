from scripts import core_simulation, plot_add_ons, test_core_simulation
import matplotlib.pyplot as plt
import numpy as np


test_core_simulation()

if __name__ == '__main__':
    # example run
    well_head_df, wetting_front_location_data = \
        core_simulation(horiz_cond=5 / 24 / 60, well_length=10, simulation_time=400, well_radius=0.2, porosity=0.3,
                        initial_VWC=0.1, recharge_duration=30, recharge_rate=30 / 60,
                        wetting_front_output_times=np.linspace(0, 400, 5))
    well_head_df.plot(x='time', y='head_in_well')
    plot_add_ons(title='Water level in the well vs. time',
                 xlabel='time [minutes]', ylabel='level [meters]')

    for wfl in wetting_front_location_data:
        df = wfl['location_df']
        plt.plot(df['R'].values, df['Z'].values, label=wfl['time'])
    plot_add_ons(title='Evolution of the wetting front',
                 xlabel='radius [meters]', ylabel='elevation above well bottom [meters]',
                 legend_title='time [min]')



