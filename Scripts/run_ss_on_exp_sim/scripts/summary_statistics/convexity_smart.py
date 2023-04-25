from scipy.spatial import ConvexHull
import numpy as np
import matplotlib.pyplot as plt
import math
import alphashape
from shapely.geometry import MultiPolygon

from Scripts.run_ss_on_exp_sim.scripts.helper_functions.helperFunctions import load_cellStates
from Scripts.run_ss_on_exp_sim.scripts.helper_functions.BacteriaFeatures import convert_um_pixel, find_bacteria_endpoints

"""
new convexity based on alpha shape, being able to choose the best alpha
"""


# alterd from BacteriaFeatures.find_bacteria_endcoordinates
def find_bacteria_endcoordinates(cs):
    """
    find the end coordinates of the cells in cs
    @param cs: dictionary containing bacteria features value

    Returns: 4 lists:   end point 1 x value
                        end point 1 y value
                        end point 2 x value
                        end point 2 y value
    """
    # end coordinates of bacteria\
    end_point1_x = [cs[it].ends[0][0] for it in cs]
    end_point1_y = [cs[it].ends[0][1] for it in cs]
    end_point2_x = [cs[it].ends[1][0] for it in cs]
    end_point2_y = [cs[it].ends[1][1] for it in cs]

    return end_point1_x, end_point1_y, end_point2_x, end_point2_y

def cal_perimeter(ordered_points):
    """
    Calculates the perimeter of some ordered list of points
    @param ordered_points numpy array
    """
    #print(type(ordered_points))
    perimeter_length  = 0

    for it in range(len(ordered_points)):
        if it == 0:
            continue
        pre = ordered_points[it-1]
        now = ordered_points[it]
        perimeter_length = perimeter_length + math.sqrt((pre[0] - now[0])**2 + (pre[1] - now[1])**2)
    perimeter_length = perimeter_length + math.sqrt((ordered_points[0][0] - ordered_points[-1][0])**2 + (ordered_points[0][1] - ordered_points[-1][1])**2)
    
    #print(perimeter_length)
    return perimeter_length


def get_boundary(cs, um_pixel_ratio=0.144, shape=0.048):
    """
    This function takes in a dictionary of cell states (cs),
    and generates a 2D NumPy array representing the filled envelope of the micro-colony.
    The envelope is filled using a morphological closing operation on a grayscale image generated from
    the input cell state dictionary.
    @param cs cellStates dict representing the cell states of a colony / micro-colony.
    @param um_pixel_ratio float representing the ratio of micrometers to pixels in the generated image.
    @param shape float shape parameter
    """
    # bacteria endpoints
    bacteria_end_point1_x, bacteria_end_point1_y, bacteria_end_point2_x, bacteria_end_point2_y = \
        find_bacteria_endpoints(cs)

    bacteria_endpoint_x = bacteria_end_point1_x + bacteria_end_point2_x
    bacteria_endpoint_y = bacteria_end_point1_y + bacteria_end_point2_y
    # convert to pixel
    bacteria_endpoint_x = convert_um_pixel(bacteria_endpoint_x, um_pixel_ratio)
    bacteria_endpoint_y = convert_um_pixel(bacteria_endpoint_y, um_pixel_ratio)

    points = [(x, y) for x, y in zip(bacteria_endpoint_x, bacteria_endpoint_y)]

    alpha_shape = alphashape.alphashape(points, shape)
    if isinstance(alpha_shape, MultiPolygon) or alpha_shape.boundary is None:
        shape = round(shape - 0.05, 3)
        alpha_shape = alphashape.alphashape(points, shape)

    coords = alpha_shape.boundary.coords.xy
    vert = [(x, y) for x, y in zip(coords[0], coords[1])]
    
    """try:
        coords = alpha_shape.boundary.coords.xy
        vert = [(x, y) for x, y in zip(coords[0], coords[1])]
    except:
        vert = []
        print(type(alpha_shape.boundary))
        for polygon in alpha_shape.geoms:
            coords = polygon.boundary.coords.xy
            vert = vert + [(x, y) for x, y in zip(coords[0], coords[1])]
            #print(vert)"""

    return vert

def polygon_area(coords):
    """
    This function calculate the area formed by a set of coordinates
    """
    n = len(coords)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += coords[i][0] * coords[j][1]
        area -= coords[j][0] * coords[i][1]
    area = abs(area) / 2.0
    return area

def cal_convexity(cs, fig_export_path, fig_name, um_pixel_ratio=0.144, shape = 0.048, display = False, smart = True):
    """
    Calculates the ratio of the convex hull perimeter and perimeter of the fitted colony boundary
    @param cs: dictionary, cell status
    @param fig_export_path: string, 
    @param fig_name: string,
    @param um_pixel_ratio: float
    @param shape: float
    @param display: bool, determine whether to display figure
    @param smart: bool, determine whether to automatically choose shape value or use the given
    """

    end_point1_x, end_point1_y, end_point2_x, end_point2_y = find_bacteria_endcoordinates(cs)
    
    um_pixel_ratio=0.144
    # convert um to pixles
    end_point1_x = np.array(end_point1_x) / um_pixel_ratio
    end_point1_y = np.array(end_point1_y) / um_pixel_ratio
    end_point2_x = np.array(end_point2_x) / um_pixel_ratio
    end_point2_y = np.array(end_point2_y) / um_pixel_ratio

    coordinates = list(zip(end_point1_x, end_point1_y))
    coordinates = coordinates + list(zip(end_point2_x, end_point2_y))
    coordinates = np.array(coordinates)
    #print(coordinates)
    hull = ConvexHull(coordinates)

    #print(hull.vertices)
    # calculated the length of the convex perimeter
    convex_hull_perimeter = cal_perimeter(coordinates[hull.vertices])

    ################### contour #####################
    # find the best alphashape value
    min_area = -1
    best_shape = -1
    best_boundary = 0
    alpha_list = []
    area_list = []
    if smart:
        for shape_i in range(1, 100, 1 ):
            shape_i = shape_i/1000
            #print(shape_i)
            alpha_list.append(shape_i)
            try:
                boundary = np.array(get_boundary(cs, um_pixel_ratio=um_pixel_ratio, shape=shape_i))
                area = polygon_area(boundary)
                area_list.append(area)
                if min_area < 0 or area < min_area:
                    min_area = area
                    best_shape = shape_i
                    best_boundary = boundary
            except Exception as e:
                area_list.append(-1)
                print("An error occurred at " + str(shape_i) + " :", e)
        #print(best_shape, " ", area)
        contour_perimeter = cal_perimeter(best_boundary)
    else:
        boundary = np.array(get_boundary(cs, um_pixel_ratio=um_pixel_ratio, shape=shape))
        contour_perimeter = cal_perimeter(boundary)


######################### plot ###########################
    # plot all endpoints of the colony
    plt.plot(coordinates[:,0], coordinates[:,1], 'o')

    #Plot convex hull:
    """for simplex in hull.simplices:
        plt.plot(coordinates[simplex, 0], coordinates[simplex, 1], 'k-')"""
    # We could also have directly used the vertices of the hull, which for 2-D are guaranteed to be in counterclockwise order:
    plt.plot(coordinates[hull.vertices,0], coordinates[hull.vertices,1], 'r--', lw=2)
    plt.plot(coordinates[hull.vertices[0],0], coordinates[hull.vertices[0],1], 'ro')
    #plt.title("convex_hull.png")
    #if display:
        #plt.show()
    #plt.close()

    # plot countour
    plt.plot(boundary[:,0], boundary[:,1], 'o')
    plt.plot(boundary[:,0], boundary[:,1], 'b--', lw=2)
    plt.title("convex_hull+boundary.png")
    if display:
        plt.show()
    plt.savefig(fig_export_path + fig_name + "_" + "hull_n_boundary.png")
    plt.close()

    # plot area respect to alpha value 
    plt.plot(alpha_list, area_list)
    if display:
        plt.show()
    plt.close()

    

    return convex_hull_perimeter / contour_perimeter, best_shape


if __name__ == '__main__':
    picklefile = "sep8_step-000089.pickle"  # "sep7_step-000097" "step-00200.pickle" #"circle_step-01000.pickle" 'jan3_step-000097.pickle' "sep8_step-000089.pickle"
    cs = load_cellStates("", picklefile)
    print(cal_convexity(cs, "", "fig_name", um_pixel_ratio=0.144, display = True))