"""
This module provides functions for image processing, including conversion between polygons and masks, dilation of masks, and extraction of polygons from masks.

Module Functions:
- polygon_to_mask(polygon, dim, style="fill", show=True): Converts a polygon to a binary mask.
- dilate_mask(mask, strength, show=True): Dilates a binary mask.
- mask_to_polygon(mask): Extracts a polygon from a binary mask.
"""

import numpy as np
from scipy.ndimage.morphology import binary_dilation
from PIL import Image, ImageDraw, ImageOps
from imantics import Mask

def polygon_to_mask(polygon, dim, style="fill", show=True):
    """
    Converts a polygon to a binary mask.

    Parameters:
    - polygon (list): List of (x, y) coordinates representing the polygon.
    - dim (tuple): Dimensions (width, height) of the output mask.
    - style (str): Style of the mask ("fill" for filled, "outline" for outlined).
    - show (bool): If True, displays the generated mask.

    Returns:
    - tuple: Binary mask as a NumPy array and corresponding PIL Image.

    Raises:
    - None
    """
    if style == "fill":
        outline, fill = None, 255
    elif style == "outline":
        outline, fill = 255, None

    img = Image.new('L', (dim[0], dim[1]), color=0)
    ImageDraw.Draw(img).polygon(polygon, outline=outline, fill=fill)
    img = ImageOps.flip(img)
    mask = np.array(img)

    if show:
        img.show()

    return mask, img

def dilate_mask(mask, strength, show=True):
    """
    Dilates a binary mask.

    Parameters:
    - mask (np.array): Binary mask as a NumPy array.
    - strength (float): Strength of dilation.
    - show (bool): If True, displays the dilated mask.

    Returns:
    - tuple: Dilated mask as a NumPy array and corresponding PIL Image.

    Raises:
    - None
    """
    strength = int(round(strength))
    strel = np.ones((strength, strength))
    dilated = binary_dilation(mask, structure=strel)
    img = Image.fromarray(dilated)

    if show:
        img.show()

    return dilated, img

def mask_to_polygon(mask):
    """
    Extracts a polygon from a binary mask.

    Parameters:
    - mask (np.array): Binary mask as a NumPy array.

    Returns:
    - list: List of (x, y) coordinates representing the extracted polygon.

    Raises:
    - None
    """
    mask = np.flip(mask, 0)
    polygons = Mask(mask).polygons()

    polygon_lengths = [len(poly) for poly in polygons.points]
    max_val = max(polygon_lengths)
    i = polygon_lengths.index(max_val)

    new_poly = polygons.points[i]

    return new_poly
