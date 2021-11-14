# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 09:57:54 2019

@author: Chris Tracy
"""

"""
This module contains functions for reading in GridFloat data from the National Elevation
Dataset (NED) and plotting the data. NOTE: As of December 2019, this file type is no longer 
an option when downloading the data. Regardless, the module is provided below.

Routine Listing
-------
read_ned_gridfloat : Read in the GridFloat files
read_header : Read in the header file associated with GridFloat files
read_flt : Read in the *.flt file that contains the elevation data
plot_elevation : Create a contour plot of the elevation data

Notes
-------
This module is intended for GridFloat files downloaded specifically from the NED archive.
You can get these files here: <https://viewer.nationalmap.gov/basic/ >.
GridFloat format specification is here: <https://www.loc.gov/preservation/digital/formats/fdd/fdd000422.shtml >.
"""

#Global packages
from types import SimpleNamespace
import numpy as np

def read_ned_gridfloat(file):
    """
    This function will read in the float and header file. 
    It will also return the lats, lons, and elevation data from
    the function.
    
    Given a filename with the appropriate path, a namespace of latitudes, longitudes, 
    and elevation data is created. The data is regularly gridded. Bad values are given NaN
    assignments. It then packs up the lats, lons, and elevation data and returns them 
    from the function. Note: the float and header file are assumed to be in the same 
    directory. The header file must also have a lower-left (ll) corner.
    
    Parameters
    ----------
    file: The file which will be read in. Both float
    and header extensions are defined from the file.
        
    Returns
    ----------
    SimpleNamespace (has several keys, some of which are optionally included):
        lats: A latitude array of nrows in the header file starting in the lower left corner
        lons: A longitude array of ncolumns in the header file starting in the lower left corner
        elevation: Flipped elevation data from the float file
    
    Example Usage
    -----------
    >>> file = 'path/to/file'
    >>> data = read_ned_gridfloat(file)
    """
    import os.path
    
    #Extract the file name and extension separately. Define header and float file paths.
    basename, ext = os.path.splitext(file)
    hdr_file = basename + '.hdr'
    flt_file = basename + '.flt'
    
    #Read in the header and float files.
    hdr = read_header(hdr_file)
    data = read_flt(flt_file)
    
    #Change any bad values in the data to NaNs.
    data[np.isclose(data, hdr['NODATA_value'])] = np.nan
    
    #Make the shape of the data a 2D array. Flip/invert the data rows.
    data.shape = (int(hdr['nrows']), int(hdr['ncols']))
    data = np.flip(data, axis=0)
    
    #Find the lat/lon ranges the data corresponds to. This
    #can be derived from some of the header attributes.
    lons = np.arange(hdr['ncols']) * hdr['cellsize'] + hdr['xllcorner']
    lats = np.arange(hdr['nrows']) * hdr['cellsize'] + hdr['yllcorner']
    
    #Pack everything up.
    return SimpleNamespace(lats = lats, lons = lons, elevation = data)

def read_header(file):
    """
    This function reads in the header file specifically.
    A dictionary with the header data is returned. First, the file
    is opened and lines are split. All of the bytes are cast as numbers.
    Then, the header values are inserted into the dictionary.
    
    Parameters
    ------------
    file: The header file to be read in.
    
    Returns
    ------------
    hdr: A dictionary with all of the desired header values.
    
    Example Usage
    ------------
    >>> file = '/path/to/file'
    >>> data = read_header(file)
    """
    
    #Use a dictionary to hold the header information.
    hdr = dict()
    
    #Open the file, iterating row-by-row.
    with open(file, 'r') as f:
        for line in f:
            entity = line.split()
            print(entity)
            
            #Most of the values are numbers, so cast them.
            #Otherwise, keep it as a string.
            if 'byteorder' in entity[0]:
                val = entity[1]
            else:
                val = float(entity[1])
            
            hdr[entity[0]] = val
        
    return hdr

def read_flt(file):
    """
    This function will read in the float file specifically.
    
    An array with the float data is returned. First, the file
    is opened and read in binary. All of the bytes are cast as numbers, as
    seen in the with open statement for the float file.
    
    Parameters
    ----------
    file: The float file to be read in.
    
    Returns
    ----------
    data: An array with all of the desired float data.
    
    Example Usage
    ----------
    >>> file = '/path/to/file'
    >>> data = read_flt(file)
    """
    
    #Read in the float file, making the data 32-bit floats.
    with open(file, 'rb') as f:
        data = np.fromfile(f, dtype = np.float32)
    
    return data

def plot_elevation(lons, lats, elevation, cmap='Reds_r', levels = None,
                   inset_bounds = None):
    """
    Given longitudes, latitudes, and elevation data, create a filled contour
    of the data. The data should be regularly gridded. A colorbar for the plot
    is also created. Optionally, you can make an inset of a region of interest.
    
    Parameters
    ----------
    lons: A 1D N-element array of longitudes of the data
    lats: A 1D M-element array of latitudes of the data
    elevation: An N x M array of the elevation data. The first
               dimension corresponds to the lons, and the second dimension
               corresponds to the lats.
    cmap: An optional string of a colormap name from Matplotlib. Default is to use 'Reds_r'.
    levels: An optional 1D array of the levels to be contoured. If None, 10 levels spanning the
            min/max of the elevation are used.
    inset_bounds: Optional list in the form of [lon_min, lon_max, lat_min, lat_max].
                  If provided, an inset plot is created.
    
    Returns
    ----------
    SimpleNamespace (Returned value has several keys, some of which are optionally included):
        contour: Matplotlib QuadContourSet
        colorbar: Matplotlib Colorbar class
        inset: Matplotlib QuadContourSet (if inset_bounds = True)
        polygon: Matplotlib Polygon (if inset_bounds = True). The polygon of the inset.
        lines: A list of Matplotlib CollectionPatches (if inset_bounds = True).
               These are the lines that connect the inset.
                        
    Example Usage
    ------------
    >>> from scipy.stats import multivariate_normal
    >>> import numpy as np
    #Generate lats and lons
    >>> lons = np.arange(100)
    >>> lats = np.arange(100)
    #Make them 2D coordinate arrays
    >>> lons2d, lats2d = np.meshgrid(lons, lats)
    #Fill in coordinates to find the 2D Gaussian PDF
    >>> coord = np.empty(lons2d.shape + (2,))
    >>> coord[:, :, 0] = lons2d
    >>> coord[:, :, 1] = lats2d
    >>> mean = [1, -2] #Coordinates of Gaussian peak
    >>> cov = [[5.0, 0], [0.0, 20.0]] #Covariance array - width of peak
    #Generate the elevation data
    >>> rv = multivariate_normal(mean, cov)
    >>> data = rv.pdf(coord)
    >>> data = data/np.max(data)*400 #Scale the data
    >>> elev_plot = plot_elevation(lons, lats, data, inset_bounds = [40, 50, 30, 40])
    
    """
    
    #Import the needed packages
    import matplotlib.pyplot as plt
    import matplotlib.colors as mplcol
    import matplotlib.ticker as mticker
    from matplotlib.patches import Polygon, ConnectionPatch
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from metpy.plots import USCOUNTIES
    
    #Define plot artist namespace
    all_plots = SimpleNamespace()
    
    #Define the default levels
    if levels is None:
        levels = np.linspace(np.nanmin(elevation), np.nanmax(elevation), 10)  
    
    #Build the colormap based on a matplotlib colormap, with two color extensions.
    #Raise an error if the input cmap cannot be found.
    try:
        colors = plt.cm.get_cmap(cmap)(np.linspace(0, 1, levels.size + 2))
    except Exception as e:
        print(e)
        return None
    
    #Set the over and under thresholds of the colormap.
    cmap_mod = mplcol.ListedColormap(colors[1:-1, :])
    cmap_mod.set_over(colors[-1, :])
    cmap_mod.set_under(colors[0, :])
    contourProps = {'levels': levels, 'cmap': cmap_mod, 'extend': 'both'}
    
    #Plot the elevation contours on the main figure. Add the contours to the namespace.
    #Add geographical features, lat/lon grid lines and axis labels.
    fig, ax = plt.subplots(figsize = (10, 6), subplot_kw = dict(projection = ccrs.PlateCarree()))
    c_big = ax.contourf(lons, lats, elevation, **contourProps)
    all_plots.contour = c_big
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(cfeature.STATES)
    ax.add_feature(cfeature.GSHHSFeature(levels = [1]))
    ax.add_feature(USCOUNTIES.with_scale('500k'), linewidth = 1.3)
    ticks = 0.25
    grid = ax.gridlines(draw_labels = True, color = 'dimgrey', linestyle = '--')
    grid.xlocator = mticker.MultipleLocator(ticks)
    grid.ylocator = mticker.MultipleLocator(ticks)
    grid.xlabels_top = False
    grid.ylabels_right = False
    
    #Make the colorbar. Add to the namespace.
    cb = fig.colorbar(c_big, ax = ax, shrink = 0.75)
    cb.set_label('Elevation (meters)', rotation = 'vertical')
    all_plots.colorbar = cb
    
    #Make an inset if inset_bounds is provided.
    if inset_bounds is not None:
        
        #Set the corners of the polygon
        lat_poly = [inset_bounds[2], inset_bounds[3]]
        lon_poly = [inset_bounds[0], inset_bounds[1]]
    
        #Using the defined corners, make the vertices for the polygon.
        lat_vert = [lat_poly[0], lat_poly[1], lat_poly[1], lat_poly[0]]
        lon_vert = [lon_poly[0], lon_poly[0], lon_poly[1], lon_poly[1]]
        lon_lat = np.column_stack((lon_vert, lat_vert))
    
        #Make the polygon patch and add to the main axis. Add to the namespace.
        poly = Polygon(lon_lat)
        c_big.ax.add_patch(poly)
        poly.set_facecolor('none')
        poly.set_edgecolor('blue')
        poly.set_linewidth(1.3)
        all_plots.polygon = poly
    
        #Define a new axis for the inset and make the inset contours.
        newax = ax.figure.add_axes([0, 0, 0.5, 0.25], projection = ccrs.PlateCarree())
        newax.outline_patch.set_edgecolor('blue')
        newax.outline_patch.set_linewidth(1.3)
        c_inset = newax.contourf(lons, lats, elevation, **contourProps)
        all_plots.inset = c_inset #Add the inset to the namespace.
        
        #Get the corners of the main axis. Set the width and height of the inset.
        ll, ul, lr, ur = ax.get_position().corners()
        inset_width = 0.35
        inset_height = 0.3
        new_pos = [lr[0]-inset_width, lr[1], inset_width, inset_height]
        c_inset.ax.set_position(new_pos)
        
        #Next, "zoom" for the inset. Modify the ticks of the inset.
        c_inset.ax.set_xlim(lon_poly)
        c_inset.ax.set_ylim(lat_poly)
        c_inset.ax.tick_params(labelbottom=False, labelleft=False)
        c_inset.ax.tick_params(bottom = False, left = False)
    
        #Add the lines. Get the lower left/upper right corners of the polygon and inset.
        #Apply the patch dictionary to the main axis for each of the lines.
        ll_poly = [lon_poly[0], lat_poly[0]]
        ur_poly = [lon_poly[1], lat_poly[1]]
        ll_inset, _, _, ur_inset = c_inset.ax.get_position().corners()
        patch_props = {'coordsA': 'data', 'axesA': c_big.ax,
                       'coordsB': 'figure fraction', 'axesB': c_inset.ax}
        line1 = ConnectionPatch(xyA=ll_poly, xyB=ll_inset, **patch_props)
        c_big.ax.add_patch(line1)
        line2 = ConnectionPatch(xyA=ur_poly, xyB=ur_inset, **patch_props)
        c_big.ax.add_patch(line2)
    
    return all_plots

# if __name__ == '__main__':
    
#     f = 'C:/Users/chris/OneDrive/Desktop/aes509/n34w112/floatn34w112_1.flt'
    
#     data = read_ned_gridfloat(f)
    
#     elv_plot = plot_elevation(data.lons, data.lats, data.elevation, inset_bounds=[-111.9, -111.4, 33.7, 33.9])