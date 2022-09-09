import sys
sys.path.append('../scripts')
from MicroColonyAnalysis import micro_colony_analysis


if __name__ == '__main__':
    input_directory = '../input files/CellProfiler_Sample'
    summary_statistic_method_list = ["Aspect Ratio", "Anisotropy", 'Density', 'dist_vs_growth_rate']
    # unit: um
    min_size_of_micro_colony = 2
    max_distance_between_cells = 3.4
    # convert micro meter to pixel
    um_pixel_ratio = 0.144
    dt = 0.05  # interval time: 3 min
    # mean_summary_statistic_report = {summary statistic name: mean value, ....}
    mean_summary_statistic_report = micro_colony_analysis(input_directory, summary_statistic_method_list, dt,
                                                     min_size_of_micro_colony, max_distance_between_cells,
                                                     um_pixel_ratio)
    # mean value of each summary statistic
    print(mean_summary_statistic_report)
