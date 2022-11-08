from scipy.spatial import ConvexHull
import imageio
import numpy as np
import matplotlib.pyplot as plt
from Scripts.SummaryStatisticsOnExpAndSim.scripts.SummaryStatistics.process import cell_data
from Scripts.SummaryStatisticsOnExpAndSim.scripts.SummaryStatistics.repo.FractalDimension.MBdimension import fractal_dimension


def save_fig(cells, hull, points, fig_export_path, fig_name):
    """
    @param cells cellStates dict
    @param hull object convex hull
    @param points numpy array bacteria end points
    @param  fig_export_path directory to export images
    @param fig_name image name

    """
    # vertexes
    vertex1_x, vertex1_y, vertex2_x, vertex2_y = cell_data.find_vertex(cells)
    minor = cell_data.find_bacteria_minor(cells)

    fig, ax = plt.subplots()
    for i in range(len(minor)):
        plt.plot([vertex1_x[i], vertex2_x[i]], [vertex1_y[i], vertex2_y[i]], lw=minor[i], solid_capstyle="round")
    for simplex in hull.simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

    fig.savefig(fig_export_path + fig_name + ".png")
    # close fig
    fig.clf()
    plt.close()

    # fills area
    fig, ax = plt.subplots()
    plt.fill(points[hull.vertices, 0], points[hull.vertices, 1], 'k')
    fig.savefig(fig_export_path + fig_name + "_fill.png")
    # close fig
    fig.clf()
    plt.close()


def create_convex_hull(cells):
    """
    goal: find convex hull
    @param cells cellStates dict
    """
    # find endpoints
    endpoints_x, endpoints_y = cell_data.bacteria_end_points(cells)

    points = np.column_stack((endpoints_x, endpoints_y))
    hull = ConvexHull(points)

    return hull, points


def calc_fractal_dimension(cells, fig_export_path, fig_name):
    """
    goal: calculation of fractal dimension (Minkowski–Bouligand dimension)

    @param  fig_export_path directory to export images
    @param fig_name image name
    @return fractal dimension
    """
    hull, points = create_convex_hull(cells)
    save_fig(cells, hull, points, fig_export_path, fig_name)
    # Import the image in greyscale
    greyscale_img = imageio.v2.imread(fig_export_path + fig_name + '_fill.png', as_gray="True") / 255.0
    return fractal_dimension(greyscale_img)
