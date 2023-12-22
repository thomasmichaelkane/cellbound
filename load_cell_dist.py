"""Load an arbitrary shaped cell distribution.

This script allows for visualisation and analysis of an arbitrary shape
cell mosaic. A distribution of cells is loaded and a CellDist object is
created for voronoi analysis.

Example
-------
Use the command line interface to load a csv containing cell coordinates
with scaling and dimensions specified in following arguments. Here we
are loading the cell distribution from the csv file dist.csv. In this
image there are 1.4 microns per pixel, and the dimensions of the image
are 6000 x 7000 pixels::

    $ python load_cell_dist.py examples/dist.csv 1.4 6000 7000

Notes
-----
    There are three required arguments denoted below (file, microns per
    pixel, and dimension (1)), but further settings with more details can
    be altered in the settings file.

Arguments
----------
file name : str
    The relative or absolute path to the coordinates file. This
    file should be a two column csv file containing x and y 
    coordinates of all cell centres (in pixels).

microns per pixel : float
    The scaling factor. How many microns per pixel in the original image.
    
dimension (1) : int
    The first (or only) dimension - X. If no other argumented are parsed,
    this singular dimension will be used to make a square. (PIXELS)
    
{optional} dimension (2) : int
    The second dimension - Y. (PIXELS)
"""

import sys

from lib.cells.cell_dist import CellDist
from lib.utils import parse
from lib.utils.utils import *

def main():
    
    FILE_NAME, MPP, DIM = parse_args()

    coords = csv_to_list(FILE_NAME)
    points = coords_as_np(coords)
    
    file_id = get_id(FILE_NAME)
    
    dist = CellDist(points, DIM, MPP, id=file_id)

    dist.show_voronoi()
    dist.find_border()
    dist.dilate_border()
    dist.print_report_full()
    dist.show_voronoi()
    
    dist.show_histograms()
    dist.save()
    dist.show_mask()
    dist.show_dil_mask()
    dist.show_points_with_border()

       
def parse_args():

    if len(sys.argv) == 1:
        
        raise KeyError("No file specified")
    
    elif len(sys.argv) == 2:
        
        raise KeyError("No scaling specified")
        
    elif len(sys.argv) == 3:
        
        raise KeyError("No dimensions specified")
        
    elif len(sys.argv) == 4:
        
        file_name = parse.file_name(sys.argv[1])
        mpp = parse.mpp(sys.argv[2])
        dim_uni = sys.argv[3]
            
        dim = [parse.dim(dim_uni), parse.dim(dim_uni)]
        
    elif len(sys.argv) == 5:
        
        file_name = parse.file_name(sys.argv[1])
        mpp = parse.mpp(sys.argv[2])
        dim = [parse.dim(sys.argv[3]), parse.dim(sys.argv[4])]   
               
    else:
        
        raise KeyError("Too many input arguments")
    
    return file_name, mpp, dim

main()