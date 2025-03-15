import logging
import numpy as np

from region_grow.representativeness_exceptions import IllegalArgumentError

class Point(object):
    """Return point (x, y) coordinates as object.
    """
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y


class RegionGrow(object):
    """Class to grow a region in an image based on given point of the matrix.
    """

    def __init__(self, image, seed_cell, threshold_shape, threshold, con):
        """Initialize the class.

        :param image: an image in the form of a np array
        :param seed_cell: the original seed list
        :param threshold_shape: shape threshold
        :param threshold: land cover threshold
        :param con: connectivity (either 4 or 8)
        """
        self.check_parameters(con, threshold_shape)
        self.image = image
        self.cell = seed_cell
        self.thresh = threshold
        self.threshold_shape = threshold_shape
        self.connectivity = con


    @staticmethod
    def check_parameters(con, threshold_shape):
        """Check the reasonability of parameters.

        :param con: connectivity
        :param threshold_shape: shape threshold
        """
        sup_cons = (4, 8)
        if con not in sup_cons:
            raise IllegalArgumentError(
                f'Connectivity value {con} not supported. Supported values '
                f'are {sup_cons}')

        shape_threshold_max = 1
        if threshold_shape > shape_threshold_max:
            raise IllegalArgumentError(
                f'Shape threshold must be {shape_threshold_max} or lower. It '
                f'is {threshold_shape} instead.'
            )

    def get_gray_dDiff(self, current_point, tmp_point):
        """Analyse if region grow is still within the same land cover.
        """
        curr_pixel = int(self.image[current_point.y, current_point.x])
        tmp_pixel = self.image[tmp_point.y, tmp_point.x]

        gray_diff = abs(int(curr_pixel) - int(tmp_pixel))

        return gray_diff

    def get_neighbours(self):
        """Get coordinates of surrounding cells based on the connectivity.

        The set of coordinates is in the direction (NW, W, SW, S, SE, E, NE, N).

        :return: tuple of coordinates in form ((x, y), (x, y), ...)
        """
        if self.connectivity == 8:
            neighbours = (
                (-1, -1), (0, -1), (1, -1), (1, 0),
                (1, 1), (0, 1), (-1, 1), (-1, 0)
            )
        else:
            neighbours = ((0, -1), (1, 0), (0, 1), (-1, 0))

        return neighbours

    def bbox2(self, img):
        """Calculate bounding box in NumPy.
        """
        rows = np.any(img, axis=1)
        cols = np.any(img, axis=0)
        rmin, rmax = np.where(rows)[0][[0, -1]]
        cmin, cmax = np.where(cols)[0][[0, -1]]
        return rmin, rmax, cmin, cmax

    def rectangularity(self, img):
        """Calculate the rectangularity measure.

        http://www.cyto.purdue.edu/cdroms/micro2/content/education/wirth10.pdf
        Fk = ratio of region area and the area of a bounding rectangle

        :param img:
        :return:
        """
        bb = self.bbox2(img)
        obj2 = np.zeros(img.shape)

        obj2[bb[0]:bb[1] + 1, bb[2]:bb[3] + 1] = 1

        # ratio of region area and the area of a bounding rectangle
        f_k = np.sum(img==1) / np.sum(obj2)

        return f_k

    def grow(self):
        """Run region grow.

        :return: binary image with the region pixels labeled
        """
        height, width = self.image.shape
        # array for accepted grown pixels
        seedmark = np.zeros((height, width))
        # LUCAS point as a seed
        seed_list = list(self.cell)

        label = 1
        label2 = 2
        # initialize seedmark with LUCAS point
        seedmark[seed_list[0].y, seed_list[0].x] = label

        neighbours_coords = self.get_neighbours()

        # rectangularity property of the object
        rect = 1.0
        while rect >= self.threshold_shape:
            next_seed_list = []
            seedmark_last = seedmark.copy()
            # new 'snake' loop around the object
            for current_point in seed_list:
                # take new point to analyze grow
                # mark the current point with the label value
                seedmark[current_point.y, current_point.x] = label

                # take a round given the connectivity selection
                for x, y in neighbours_coords:
                    x_tmp = current_point.x + x
                    y_tmp = current_point.y + y

                    # check position within tile patch
                    if x_tmp < 0 or y_tmp < 0 or x_tmp >= width or y_tmp >= height:
                        continue

                    # checking if the RG can grow
                    gray_diff = self.get_gray_dDiff(
                        current_point, Point(x_tmp, y_tmp)
                    )

                    if gray_diff < self.thresh and seedmark[y_tmp, x_tmp] == 0:
                        # if new valid point, append it
                        seedmark[y_tmp, x_tmp] = label
                        new_point = Point(x_tmp, y_tmp)

                        next_seed_list.append(new_point)

                    elif gray_diff >= self.thresh and seedmark[y_tmp, x_tmp] == 0:
                        # edge to other class
                        seedmark[y_tmp, x_tmp] = label2

            rect_last = rect
            rect = self.rectangularity(seedmark)

            # if no more cells to go, stop
            if len(next_seed_list) == 0:
                break

            seed_list = next_seed_list

        if rect_last > 1.0:
            rect_last = 1.0

        logging.debug("End of GROW!")
        logging.debug(f'Rectangularity: {rect_last}')
        logging.debug(f"Size of object: {np.sum(seedmark_last)} (px).")

        return seedmark_last, rect_last
