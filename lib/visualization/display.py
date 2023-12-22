"""
This module provides functions for visualizing data, including histograms, point plots, Voronoi diagrams, and alpha plots.

Module Functions:
- histogram(pop, xlabel, bins=15): Displays a histogram of a given population.
- plot_text(x, y, text, color, size): Plots text at a specified position on the graph.
- plot_point(x, y, index=None, neighbors=None, style='.', color="r", size="3"): Plots a point on the graph.
- plot_sides(polygon): Plots the sides of a polygon on the graph.
- fill_cell(polygon, color): Fills the interior of a polygon with a specified color.
- show_voronoi(vor, cells): Displays a Voronoi diagram.
- add_neighbor_legend(ax): Adds a legend for the number of neighbors in a Voronoi diagram.
- get_neighbors_color(num_neighbors): Calculates the color based on the number of neighbors.
- show_polygon(polygon_vertices, dim, points=None): Displays a polygon on the graph.
- alpha_plots(edge_sets, points, alphas): Displays alpha plots based on edge sets and points.
"""

import matplotlib.pyplot as plt
from matplotlib import cm, colors, patches

from ..utils.utils import *
from ..utils.enums import PointType
from ..utils.settings import display_settings

def histogram(pop, xlabel, bins=15):
    """
    Displays a histogram of a given population.

    Parameters:
    - pop (dict): Population data in dictionary format.
    - xlabel (str): Label for the x-axis.
    - bins (int): Number of bins in the histogram.

    Returns:
    - None
    """
    pop_list = list(pop.values())
    plt.hist(pop_list, bins=bins)
    plt.xlabel(xlabel)
    plt.show()

def plot_text(x, y, text, color, size):
    """
    Plots text at a specified position on the graph.

    Parameters:
    - x (float): X-coordinate of the text position.
    - y (float): Y-coordinate of the text position.
    - text (str): Text to be displayed.
    - color (str): Color of the text.
    - size (int): Font size of the text.

    Returns:
    - None
    """
    plt.text(x, y, text, fontdict=None, ha="center", va="center", color=color, fontsize=size)

def plot_point(x, y, index=None, neighbors=None, type=PointType.NONE, style='.', color="r", size="3"):
    """
    Plots a point on the graph.

    Parameters:
    - x (float): X-coordinate of the point.
    - y (float): Y-coordinate of the point.
    - index (int): Index value for the point.
    - neighbors (int): Number of neighbors for the point.
    - type (PointType): Type of point.
    - style (str): Marker style for the point.
    - color (str): Color of the point.
    - size (int): Size of the point.

    Returns:
    - None
    """
    if type == PointType.INDEX:
        plot_text(x, y, index, color, size)
    elif type == PointType.NEIGHBORS:
        plot_text(x, y, neighbors, color, size)
    else:
        plt.plot(x, y, marker=style, color=color, markersize=size)

def plot_sides(polygon):
    """
    Plots the sides of a polygon on the graph.

    Parameters:
    - polygon (list): List of (x, y) coordinates representing the polygon.

    Returns:
    - None
    """
    first_vertex = polygon[0]
    polygon.append(first_vertex)
    plt.plot(*zip(*polygon), color=display_settings["line_color"], 
                             linewidth=display_settings["line_width"])

def fill_cell(polygon, color):
    """
    Fills the interior of a polygon with a specified color.

    Parameters:
    - polygon (list): List of (x, y) coordinates representing the polygon.
    - color (str): Color for filling the polygon.

    Returns:
    - None
    """
    plt.fill(*zip(*polygon), color)

def show_voronoi(vor, cells):
    """
    Displays a Voronoi diagram.

    Parameters:
    - vor: Voronoi diagram object.
    - cells (list): List of indices representing Voronoi cells.

    Returns:
    - None
    """
    ax = plt.axes()

    for j in cells:
        region = vor.regions[vor.point_region[j]]
        (x, y) = vor.points[j]
        cell_polygon = [vor.vertices[i] for i in region]
        num_neighbors = len(cell_polygon)

        if display_settings["line_style"] is not None:
            plot_sides(cell_polygon)

        if display_settings["fill"]:
            color = get_neighbors_color(num_neighbors)
            fill_cell(cell_polygon, color)

        if display_settings["point_type"] is not PointType.NONE:
            plot_point(x, y, index=j,
                       neighbors=num_neighbors,
                       type=display_settings["point_type"],
                       style=display_settings["point_style"],
                       color=display_settings["point_color"],
                       size=display_settings["point_size"])

    if display_settings["fill"]:
        add_neighbor_legend(ax)

    plt.show()

def add_neighbor_legend(ax):
    """
    Adds a legend for the number of neighbors in a Voronoi diagram.

    Parameters:
    - ax: Matplotlib axis object.

    Returns:
    - None
    """
    legend_patches = []
    neighbor_range = inc_range(*display_settings["scale_neighbors"])

    for i in neighbor_range:
        color = get_neighbors_color(i)
        patch = patches.Patch(color=color, label=i)
        legend_patches.append(patch)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax.legend(handles=legend_patches, loc="right", bbox_to_anchor=(1.2, 0.5), fancybox=True, shadow=True)

def get_neighbors_color(num_neighbors):
    """
    Calculates the color based on the number of neighbors.

    Parameters:
    - num_neighbors (int): Number of neighbors for a Voronoi cell.

    Returns:
    - str: Hex color code.
    """
    scale_neighbors = display_settings["scale_neighbors"]
    alpha = (num_neighbors - scale_neighbors[0]) * (1 / (scale_neighbors[1] - scale_neighbors[0]))
    cmap = cm.get_cmap(display_settings["neighbor_colormap"])
    color = colors.to_hex(cmap(alpha))
    return color

def show_polygon(polygon_vertices, dim, points=None):
    """
    Displays a polygon on the graph.

    Parameters:
    - polygon_vertices (list): List of (x, y) coordinates representing the polygon.
    - dim (tuple): Dimensions of the graph.
    - points (list): List of (x, y) coordinates representing points to be plotted.

    Returns:
    - None
    """
    fig, ax = plt.subplots()

    polygon = patches.Polygon(polygon_vertices,
                              closed=True,
                              fill=True,
                              color=display_settings["mask_color"])

    if points is not None:
        for (x, y) in points:
            plot_point(x, y,
                       color=display_settings["mask_point_color"],
                       size=display_settings["mask_point_size"])

    ax.add_patch(polygon)
    ax.set_xlim([0, dim[0]])
    ax.set_ylim([0, dim[1]])
    plt.show()

def alpha_plots(edge_sets, points, alphas):
    """
    Displays alpha plots based on edge sets and points.

    Parameters:
    - edge_sets (list): List of sets representing edges.
    - points (np.array): Array of shape (n, 2) representing points.
    - alphas (list): List of alpha values.

    Returns:
    - None
    """
    fig, ax = plt.subplots(figsize=(18, 4))
    num_plots = len(alphas)

    for k, alpha in enumerate(alphas):
        plt.subplot(1, num_plots, k + 1)
        plt.plot(points[:, 0], points[:, 1], '.')
        for i, j in edge_sets[k]:
            plt.plot(points[[i, j], 0], points[[i, j], 1])

        plt.text(2000, 6200, f"alpha={alpha}", size=18)

    plt.show()
