#from Scripts.RunSS_henry.scripts.process.BacteriaFeatures import fit_ellipse
import numpy as np
import pandas as pd

from Scripts.run_ss_on_exp_sim.scripts.lownerJohnEllipseMaster.src.lownerJohnEllipse import welzl


def fit_enclosing_ellipse(points):
    """
    goal: An ellipse is fitted around bacteria in a specific micro colony (local mode) or
    in specific time step (global mode).

    @param points dataframe endpoints (vertex & co-vertex) of bacteria in specific micro colony(local mode) or
    in specific time step (global mode)

    Returns: enclosing_ellipse tuple a tuple (c, a, b, t), where c = (x, y) is the center, a and
            b are the major and minor radii, and t is the rotation angle.

    """

    # convert dataframe to numpy array
    points = points.to_numpy()
    # print(points)
    # finds the smallest ellipse covering a finite set of points
    # https://github.com/dorshaviv/lowner-john-ellipse
    enclosing_ellipse = welzl(points)

    return enclosing_ellipse

def fit_ellipse(cs):
    """

    @param cs dict bacteria features value in specific micro colony(local mode) or
    in specific time step (global mode)

    Returns: ellipse_params  tuple a tuple (c, a, b, t), where c = (x, y) is the center, a and
            b are the major and minor radii, and t is the rotation angle.

    """
    bacteria_center_x = [cs[it].pos[0] for it in cs]
    bacteria_center_y = [cs[it].pos[1] for it in cs]

    endpoints = pd.DataFrame(zip(bacteria_center_x, bacteria_center_y))

    # fit ellipse
    ellipse_params = fit_enclosing_ellipse(endpoints)

    return ellipse_params

def calc_aspect_ratio(cs):
    """
    goal: calculation of Aspect Ratio
    @param cs dict bacteria features value in specific micro colony(local mode) or in specific time step (global mode)

    Returns: aspect_ratio float aspect ratio of micro colony(local mode) or specific time step (global mode)

    """

    # fit ellipse
    # if  number of bacteria < 2 : ellipse_params = `nan`
    if len(cs) > 2:
        ellipse_params = fit_ellipse(cs)

        # calculate aspect ratio
        center_pos, major, minor, theta = ellipse_params
        aspect_ratio = minor / major
    else:
        aspect_ratio = np.nan

    return aspect_ratio

def calc_aspect_ratio(cs):
    """
    goal: calculation of Aspect Ratio
    @param ellipse_params  tuple a tuple (c, a, b, t), where c = (x, y) is the center, a and
            b are the major and minor radii, and t is the rotation angle.

    Returns: aspect_ratio float aspect ratio of micro colony(local mode) or specific time step (global mode)

    """
    ellipse_params = fit_ellipse(cs)
    # calculate aspect ratio
    center_pos, major, minor, theta = ellipse_params
    aspect_ratio = minor / major

    return aspect_ratio
