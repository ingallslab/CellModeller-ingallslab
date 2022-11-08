from Scripts.SummaryStatisticsOnExpAndSim.scripts.SummaryStatistics.anisotropy_stat import anisotropy_calc
from Scripts.SummaryStatisticsOnExpAndSim.scripts.SummaryStatistics.aspect_ratio_stat import aspect_ratio_calc
from Scripts.SummaryStatisticsOnExpAndSim.scripts.SummaryStatistics import density_calculation, growth_rate_stats
from Scripts.SummaryStatisticsOnExpAndSim.scripts.SummaryStatistics.fractal_dimension import calc_fractal_dimension


def calc_summary_statistics(cs, stats_value, summary_statistic_method, max_distance_between_cells, dt, fig_export_path,
                            fig_name):
    """
    goal: call related functions to calculate summary statistics
    @param cs dict cellStates dict
    @param stats_value dict store summary statistics over time-lapse (for each individual time-step)
    @param summary_statistic_method   str     the method that we apply over time-lapse
    @param max_distance_between_cells float There is a maximum distance between bacteria that can be neighbours.
    @param dt float interval time
    @param fig_export_path str directory to export images
    @param fig_name str  image name
    @return stats_value dict summary statistics values over time-lapse (for each individual time-step)
    """

    # calculation of aspect ratio
    if "Aspect Ratio" in summary_statistic_method:
        aspect_ratio = aspect_ratio_calc(cs)
        # store aspect ratio
        stats_value["Aspect Ratio"].append(aspect_ratio)

    # calculation of Anisotropy
    if "Anisotropy" in summary_statistic_method:
        mean_anisotropy = anisotropy_calc(cs, max_distance_between_cells)
        # store anisotropy
        stats_value["Anisotropy"].append(mean_anisotropy)

    # calculation of Density
    if "Density" in summary_statistic_method:
        density = density_calculation.main(cs)
        # store anisotropy
        stats_value["Density"].append(density)

    # calculation of Correlate growth penalty based on location
    # Note that this summary stat currently is calc. only for time steps containing more than 2 bacteria
    # because the center of the bacteria is used to fit ellipse but the ellipse can be fitted when at
    # least three points are available
    if "dist_vs_growth_rate" in summary_statistic_method and len(cs) > 2:
        dist_vs_growth_rate = growth_rate_stats.main(cs, dt)
        # store dist_vs_growth_rate
        stats_value["dist_vs_growth_rate"].append(dist_vs_growth_rate)

    # calculation of fractal dimension
    if "fractal dimension" in summary_statistic_method:
        fractal_dimension = calc_fractal_dimension(cs, fig_export_path, fig_name)
        # store fractal dimension
        stats_value["fractal dimension"].append(fractal_dimension)

    return stats_value
