import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import AMRh5
import matplotlib.animation as animation
# import xarray as xr
# import cartopy.crs as ccrs

def plot_map(flatAMRh5obj, kwargs_map={}, kwargs_colorbar={}, **kwargs_pcolormesh):
    """
    Plots a 2D map from an AMRh5 object using matplotlib.
    
    Parameters:
    - AMRh5_obj: An instance of the AMRh5 class containing the data to plot.
    - kwargs_map: Dictionary of keyword arguments for map customization (e.g., title, labels).
    - kwargs_colorbar: Dictionary of keyword arguments for colorbar customization.
    - kwargs_imshow: Additional keyword arguments for plt.imshow().
    
    Returns:
    - fig, ax: The matplotlib figure and axis objects.
    """
    
    x = flatAMRh5obj.x
    y = flatAMRh5obj.y
    data = flatAMRh5obj.data
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot the data using pcolormesh
    cax = ax.pcolormesh(x, y, data, shading='auto', **kwargs_pcolormesh)
    
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


def plot_diff_map(flatAMRh5obj1, flatAMRh5obj2, kwargs_map={}, kwargs_colorbar={}, **kwargs_pcolormesh):

    """
    Plots the difference between two 2D maps from AMRh5 objects using matplotlib.
    
    Parameters:
    - AMRh5_obj1: First instance of the AMRh5 class.
    - AMRh5_obj2: Second instance of the AMRh5 class.
    - kwargs_map: Dictionary of keyword arguments for map customization (e.g., title, labels).
    - kwargs_colorbar: Dictionary of keyword arguments for colorbar customization.
    - kwargs_imshow: Additional keyword arguments for plt.imshow().
    
    Returns:
    - fig, ax: The matplotlib figure and axis objects.
    """
    
    x = flatAMRh5obj1.x
    y = flatAMRh5obj1.y
    data_diff = flatAMRh5obj2.data - flatAMRh5obj1.data
    
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
