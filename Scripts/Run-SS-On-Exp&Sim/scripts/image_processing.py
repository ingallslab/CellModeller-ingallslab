import cv2
import numpy as np
from PIL import Image, ImageDraw
import skimage.morphology as morphology
import BacteriaFeatures


def draw_gray_scale_img(img_dimensions, bacteria_center_x, bacteria_center_y, bacteria_major, bacteria_radius,
                        bacteria_orientation):
    """
    gaol: draw gray scale bacteria image

    @param  img_dimensions   list    image dimension [width, height] (pixel)
    @param  bacteria_center_x  list    bacteria x coordinate of center position (pixel)
    @param  bacteria_center_y  list    bacteria y coordinate of center position (pixel)
    @param  bacteria_major     list    bacteria major axis length (pixel)
    @param  bacteria_radius         list    bacteria minor axis length (pixel)
    @param  bacteria_orientation   list    bacteria orientation (radian)
    @return gray_scale_img                 gray scale bacteria image
    """

    # L mode: black and white
    gray_scale_img = Image.new(mode='L', size=img_dimensions)
    draw = ImageDraw.Draw(gray_scale_img)

    # Adjust orientations to image coordinate system
    img_orientation = -np.array(bacteria_orientation)
    img_orientation_norm = img_orientation + np.pi / 2

    # lengths
    half_length_along_axis_x = np.array(bacteria_major) / 2 * np.cos(img_orientation)
    half_length_along_axis_y = np.array(bacteria_major) / 2 * np.sin(img_orientation)
    radius_perpendicular_to_axis_x = np.array(bacteria_radius) * np.cos(img_orientation_norm)
    radius_perpendicular_to_axis_y = np.array(bacteria_radius) * np.sin(img_orientation_norm)

    # convert to numpy array
    bacteria_center_x = np.array(bacteria_center_x)
    bacteria_center_y = np.array(bacteria_center_y)

    # vertices
    p1_p4_x_second_part = half_length_along_axis_x + radius_perpendicular_to_axis_x
    p1_p4_y_second_part = half_length_along_axis_y + radius_perpendicular_to_axis_y
    p2_p3_x_second_part = half_length_along_axis_x - radius_perpendicular_to_axis_x
    p2_p3_y_second_part = half_length_along_axis_y - radius_perpendicular_to_axis_y

    p5_p7_x_second_part = half_length_along_axis_x - bacteria_radius
    p5_p7_y_second_part = half_length_along_axis_y - bacteria_radius
    p6_p8_x_second_part = half_length_along_axis_x + bacteria_radius
    p6_p8_y_second_part = half_length_along_axis_y + bacteria_radius

    # rectangle
    p1 = (bacteria_center_x + p1_p4_x_second_part, bacteria_center_y + p1_p4_y_second_part)
    p2 = (bacteria_center_x + p2_p3_x_second_part, bacteria_center_y + p2_p3_y_second_part)
    p3 = (bacteria_center_x - p2_p3_x_second_part, bacteria_center_y - p2_p3_y_second_part)
    p4 = (bacteria_center_x - p1_p4_x_second_part, bacteria_center_y - p1_p4_y_second_part)

    # bacteria ends
    p5 = (bacteria_center_x + p5_p7_x_second_part, bacteria_center_y + p5_p7_y_second_part)
    p6 = (bacteria_center_x + p6_p8_x_second_part, bacteria_center_y + p6_p8_y_second_part)
    p7 = (bacteria_center_x - p5_p7_x_second_part, bacteria_center_y - p5_p7_y_second_part)
    p8 = (bacteria_center_x - p6_p8_x_second_part, bacteria_center_y - p6_p8_y_second_part)

    # start & end angle should be in degree
    end_angle = BacteriaFeatures.convert_radian_degree(img_orientation_norm)
    start_angle = end_angle - 180

    for index, value in enumerate(p1[0]):
        draw.polygon([(value, p1[1][index]), (p2[0][index], p2[1][index]), (p3[0][index], p3[1][index]),
                      (p4[0][index], p4[1][index])], fill=255, outline=None)
        draw.pieslice([(p5[0][index], p5[1][index]), (p6[0][index], p6[1][index])], start=start_angle[index],
                      end=end_angle[index], fill=255, outline=None)
        draw.pieslice([(p7[0][index], p7[1][index]), (p8[0][index], p8[1][index])], start=end_angle[index],
                      end=start_angle[index], fill=255, outline=None)

    return gray_scale_img


def save_fig(contours, fig_export_path, fig_name):

    contours_from_array = Image.fromarray(contours)
    contours_from_array.save(fig_export_path + fig_name + ".png")


def fill_contours(cs, img_dimension=None, um_pixel_ratio=0.144, margin=10):
    """
    @param cs cellStates dict
    @param margin int pixel
    """
    # bacteria features
    bacteria_center_x, bacteria_center_y = BacteriaFeatures.find_bacteria_center_position(cs)
    bacteria_center_x = BacteriaFeatures.convert_um_pixel(bacteria_center_x, um_pixel_ratio)
    bacteria_center_y = BacteriaFeatures.convert_um_pixel(bacteria_center_y, um_pixel_ratio)

    if min(bacteria_center_x) < 0 or min(bacteria_center_y) < 0:
        bacteria_center_x += (abs(min(bacteria_center_x)) + margin)
        bacteria_center_y += (abs(min(bacteria_center_y)) + margin)

    if img_dimension is not None:
        bacteria_center_x += (img_dimension[0] / 2 - np.average(bacteria_center_x))
        bacteria_center_y += (img_dimension[1] / 2 - np.average(bacteria_center_y))
    else:
        img_dimension = [int(np.ceil(1.5 * max(bacteria_center_x))), int(np.ceil(1.5 * max(bacteria_center_y)))]

        if min(bacteria_center_x) - margin > margin and min(bacteria_center_y) - margin > margin:
            bacteria_center_x -= margin
            bacteria_center_y -= margin
        else:
            bacteria_center_x += margin
            bacteria_center_y += margin


    bacteria_orientation = BacteriaFeatures.find_bacteria_orientation(cs)
    bacteria_major = BacteriaFeatures.convert_um_pixel(BacteriaFeatures.find_bacteria_major(cs))
    bacteria_minor = BacteriaFeatures.convert_um_pixel(BacteriaFeatures.find_bacteria_minor(cs))

    # Draw black and white image
    gray_scale_img = draw_gray_scale_img(img_dimension, bacteria_center_x, bacteria_center_y, bacteria_major, bacteria_minor,
                                 bacteria_orientation)

    # Perform morphological closing to get rid of any gaps
    img_cloe = morphology.closing(gray_scale_img, footprint=np.ones((7, 7), dtype='uint8'))

    # Get contours
    contours, hierarchy = cv2.findContours(np.copy(img_cloe), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Fill contours so that the image is a filled envelope of the micro-colony
    filled_img = cv2.drawContours(img_cloe, contours, -1, 255, cv2.FILLED)

    return filled_img


def find_external_contours(img, fig_export_path, fig_name):

    # Read in the image as grayscale
    im = cv2.imread(img, 0)

    # we want to find the best contour possible with CHAIN_APPROX_NONE
    contours, hierarchy = cv2.findContours(im.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

    # Create an output of all zeroes that has the same shape as the input
    boundary = np.zeros_like(im)

    cv2.drawContours(boundary, contours, -1, 255, 1)

    if fig_export_path:
        save_fig(boundary, fig_export_path, fig_name)

    return boundary


def pixel_to_coordinate(contours):

    x_coordinates = []
    y_coordinates = []
    for x_index, contour in enumerate(contours):
        white_pixel = np.where(contour == 255)
        if white_pixel[0].size > 0:
            x_coordinates.extend([x_index] * white_pixel[0].size)
            y_coordinates.extend(list(white_pixel[0]))

    return x_coordinates, y_coordinates


def coordinate_sampling(x_coordinates, y_coordinates, skip):

    x_coordinates = [x for i, x in enumerate(x_coordinates) if i % skip == 0]
    y_coordinates = [y for i, y in enumerate(y_coordinates) if i % skip == 0]
    return x_coordinates, y_coordinates

