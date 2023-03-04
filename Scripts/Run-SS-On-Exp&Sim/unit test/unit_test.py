import sys
sys.path.append('../scripts')
from DataAnalysis import data_analysis

# Note: this is a single code that is capable of calculating summary statistics of both exp. and sim.

if __name__ == '__main__':
    # cellprofiler output pickle files
    # input_directory = '../input files/CellProfiler_Sample'
    #  cellModeller output pickle files
    input_directory = '../input files/CellModeller_Sample'

    summary_statistic_method_list = ["Aspect Ratio", "Anisotropy", "fourier_descriptor"]
    # Calculation of summary stats on only the last time step
    only_last_time_step = True
    dt = 0.05  # interval time: 3 min
    # useful for calculation of fractal dimension (optional)
    fig_export_path = None

    # optional arguments
    # 1. mode: local (calculate for each micro colony) or global(calculate over time lapse without considering
    # micro-colonies) (default: 'global')
    mode = 'global'

    # 2. minimum size of micro colony (unit: um - default: 2)
    min_size_of_micro_colony = 2

    # 3. maximum distance between bacteria that can be neighbours. (unit: pixel - default: 3.4)
    max_distance_between_cells = 3.4

    # 4. convert micro meter to pixel (default: 0.144)
    um_pixel_ratio = 0.144

    # result: mean_summary_statistic_report = {summary statistic name: mean value, ....}
    # local summary statistics
    # local_mean_summary_statistic_report = data_analysis(input_directory, summary_statistic_method_list, dt,
    #                                                    mode, max_distance_between_cells, um_pixel_ratio,
    #                                                    min_size_of_micro_colony, fig_path=fig_export_path,
    #                                                          only_last_time_step=only_last_time_step)
    # global summary statistics
    global_mean_summary_statistic_report = data_analysis(input_directory, summary_statistic_method_list, dt, mode,
                                                         max_distance_between_cells, um_pixel_ratio,
                                                         fig_path=fig_export_path,
                                                         only_last_time_step=only_last_time_step)
    # mean value of each summary statistic
    # print('local mode:')
    # print(local_mean_summary_statistic_report)
    print('global mode:')
    print(global_mean_summary_statistic_report)
