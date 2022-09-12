import sys

sys.path.append('../scripts')
from DataAnalysis import data_analysis

# Note: this is a single code that is capable of calculating summary statistics of both exp. and sim.

'''
if __name__ == '__main__':
    #  cellModeller output pickle files
    input_directory = '../input files/CellModeller_Sample'
    summary_statistic_method_list = ["Aspect Ratio", "Anisotropy", 'Density', 'dist_vs_growth_rate']
    # mode: local (calculate for each micro colony) or global(calculate over time lapse without considering 
    # micro-colonies)
    mode1 = "local"
    mode2 = 'global'
    # unit: um
    min_size_of_micro_colony = 2
    max_distance_between_cells = 3.4
    # convert micro meter to pixel
    um_pixel_ratio = 0.144
    dt = 0.05  # interval time: 3 min
    # mean_summary_statistic_report = {summary statistic name: mean value, ....}
    # local summary statistics
    local_mean_summary_statistic_report = data_analysis(input_directory, summary_statistic_method_list, mode1,
                                                        dt, max_distance_between_cells, um_pixel_ratio,
                                                        min_size_of_micro_colony)
    # global summary statistics
    global_mean_summary_statistic_report = data_analysis(input_directory, summary_statistic_method_list, mode2,
                                                         dt, max_distance_between_cells, um_pixel_ratio)
    # mean value of each summary statistic
    print(local_mean_summary_statistic_report)
    print(global_mean_summary_statistic_report)
'''

if __name__ == '__main__':
    # cellprofiler output pickle files
    input_directory = '../input files/CellProfiler_Sample'
    summary_statistic_method_list = ["Aspect Ratio", "Anisotropy", 'Density', 'dist_vs_growth_rate']
    # mode: local (calculate for each micro colony) or global(calculate over time lapse without considering
    # micro-colonies)
    mode1 = "local"
    mode2 = 'global'
    # unit: um
    min_size_of_micro_colony = 2
    max_distance_between_cells = 3.4
    # convert micro meter to pixel
    um_pixel_ratio = 0.144
    dt = 0.05  # interval time: 3 min
    # mean_summary_statistic_report = {summary statistic name: mean value, ....}
    # local summary statistics
    local_mean_summary_statistic_report = data_analysis(input_directory, summary_statistic_method_list, mode1,
                                                        dt, max_distance_between_cells, um_pixel_ratio,
                                                        min_size_of_micro_colony)
    # global summary statistics
    global_mean_summary_statistic_report = data_analysis(input_directory, summary_statistic_method_list, mode2,
                                                         dt, max_distance_between_cells, um_pixel_ratio)
    # mean value of each summary statistic
    print(local_mean_summary_statistic_report)
    print(global_mean_summary_statistic_report)
