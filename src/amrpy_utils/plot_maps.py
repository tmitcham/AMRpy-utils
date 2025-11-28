import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.animation as animation

import AMRh5
from constants import OCEAN_DENSITY, ICE_DENSITY

# import xarray as xr
# import cartopy.crs as ccrs

def plot_map(flatAMRh5obj, kwargs_map=None, kwargs_colorbar=None, kwargs_pcolormesh=None):
    """
    Plots a 2D map from an AMRh5 object using matplotlib.
    
    Parameters:
    - flatAMRh5obj:      An instance of the AMRh5 class containing the data to plot.
    - kwargs_map:        Dictionary of keyword arguments for map customization (e.g., title, labels).
    - kwargs_colorbar:   Dictionary of keyword arguments for colorbar customization.
    - kwargs_pcolormesh: Additional keyword arguments for plt.pcolormesh().
    
    Returns:
    - fig, ax: The matplotlib figure and axis objects.
    """
    
    x = flatAMRh5obj.x
    y = flatAMRh5obj.y
    data = flatAMRh5obj.data
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot the data using pcolormesh
    cax = ax.pcolormesh(x, y, data, shading='auto', **(kwargs_pcolormesh or {}))
    
    # Add colorbar
    cbar = fig.colorbar(cax, ax=ax, **(kwargs_colorbar or {}))
    
    # Customize map with provided arguments
    if kwargs_map is not None:
        if 'title' in kwargs_map:
            ax.set_title(kwargs_map['title'])
        if 'xlabel' in kwargs_map:
            ax.set_xlabel(kwargs_map['xlabel'])
        if 'ylabel' in kwargs_map:
            ax.set_ylabel(kwargs_map['ylabel'])
        if 'xlim' in kwargs_map:
            ax.set_xlim(kwargs_map['xlim'])
        if 'ylim' in kwargs_map:
            ax.set_ylim(kwargs_map['ylim'])
    
    return fig, ax


def plot_diff_map(flatAMRh5obj1, flatAMRh5obj2, kwargs_map=None, kwargs_colorbar=None, kwargs_pcolormesh=None):

    """
    Plots the difference between two 2D maps from AMRh5 objects using matplotlib.

    Parameters:
    - AMRh5_obj1:        First instance of the AMRh5 class.
    - AMRh5_obj2:        Second instance of the AMRh5 class.
    - kwargs_map:        Dictionary of keyword arguments for map customization (e.g., title, labels).
    - kwargs_colorbar:   Dictionary of keyword arguments for colorbar customization.
    - kwargs_pcolormesh: Additional keyword arguments for plt.pcolormesh().
    
    Returns:
    - fig, ax: The matplotlib figure and axis objects.
    """
    
    x = flatAMRh5obj1.x
    y = flatAMRh5obj1.y
    data_diff = flatAMRh5obj2.data - flatAMRh5obj1.data

    b
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot the difference data using pcolormesh
    cax = ax.pcolormesh(x, y, data_diff, shading='auto', **kwargs_pcolormesh)
    
    # Add colorbar
    cbar = fig.colorbar(cax, ax=ax, **kwargs_colorbar)
    
    # Customize map with provided arguments
    if 'title' in kwargs_map:
        ax.set_title(kwargs_map['title'])
    if 'xlabel' in kwargs_map:
        ax.set_xlabel(kwargs_map['xlabel'])
    if 'ylabel' in kwargs_map:
        ax.set_ylabel(kwargs_map['ylabel'])
    if 'xlim' in kwargs_map:
        ax.set_xlim(kwargs_map['xlim'])
    if 'ylim' in kwargs_map:
        ax.set_ylim(kwargs_map['ylim'])
    
    return fig, ax


def plot_gls(fig, ax, glAMR5objs, kwargs_contour=None):
    """
    Plots grounding line contours on an existing map.
    
    Parameters:
    - fig:             The matplotlib figure object.
    - ax:              The matplotlib axis object.
    - glAMR5objs:      A list of AMRh5 objects containing the bedrock data, and then 
                       the thickness dataset(s) for the grounding line(s) to be plotted.
    - kwargs_contour:  Additional keyword arguments for plt.contour().
    
    Returns:
    - contour: The contour object created.
    """
    
    x = glAMR5objs[0].x
    y = glAMR5objs[0].y
    b = glAMR5objs[0].data
    h = [None] * len(glAMR5objs[1:])
    haf = h.copy()
    gl = h.copy()

    if kwargs_contour is None:
        kwargs_contour = {}


    for i in range(len(glAMR5objs[1:])):

        h[i] = glAMR5objs[i+1].data
        h[i][h[i]==0] = np.nan
        haf[i] = (-1*b)*(OCEAN_DENSITY/ICE_DENSITY)
        gl[i] = h[i] - haf[i]

        color_i = kwargs_contour.get('colors', [None])[i] if 'colors' in kwargs_contour else 'black'
        linewidth_i = kwargs_contour.get('linewidths', np.ones(len(glAMR5objs[1:])))[i] if 'linewidths' in kwargs_contour else 1.

        contour = ax.contour(x, y, gl[i], levels=[0], colors=color_i, linewidths=linewidth_i)
       
    return contour
