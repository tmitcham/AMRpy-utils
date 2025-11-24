import os
import pickle
import cv2
import matplotlib.colors as col
from matplotlib import animation
from matplotlib import pyplot as plt

import AMRh5
from plot_maps import plot_map
import create_video


# load BISICLES h5 file
print("Loading BISICLES h5 file for thickness variable...")
bisicles_file_dHdt = AMRh5.BISICLESh5("/mnt/c/Users/tm17544/OneDrive - University of Bristol/Projects/UKESM/UKESM1p3_initial_states/modern_AIS/hist_forcing/plot.modern.hist.000905.2d.hdf5", 'dThickness/dt')
flat_dHdt = bisicles_file_dHdt.flatten(lev=-1)

print("Plotting dH/dt map...")
# plot dH/dt map
fig, ax = plot_map(flat_dHdt,
                   kwargs_map={'title':f"dH/dt in year: {flat_dHdt.time}",
                               'xlabel':'X (m)', 'ylabel':'Y (m)'},
                   kwargs_colorbar={'label':'dH/dt (m/a)'},
                   cmap='RdBu',
                   norm= col.Normalize(-2, 2))
fig.savefig('dH_dt_map.png')


# Plot every file in a directory and then create a video from the figure objects

FIGS_PLOTTED = False

if not FIGS_PLOTTED:

    print("Plotting dH/dt evolution from all files in directory...")
    data_for_anim = []

    data_directory = "/mnt/c/Users/tm17544/OneDrive - University of Bristol/Projects/UKESM/UKESM1p3_initial_states/modern_AIS/hist_forcing"

    for filename in sorted(os.listdir(data_directory)):
        if filename.endswith('.hdf5') and filename.startswith('plot.'):
            print(f"Processing file: {filename}")
            
            bisicles_data= AMRh5.BISICLESh5(os.path.join(data_directory, filename), 'dThickness/dt').flatten(lev=-1)
            
            print("Plotting thickness map...")
            # plot thickness map
            fig, ax = plot_map(bisicles_data, 
                            kwargs_map={'title':f"dH/dt in year: {bisicles_data.time}",
                                        'xlabel':'X (m)', 'ylabel':'Y (m)'},
                            kwargs_colorbar={'label':'dH/dt (m/a)'}, 
                            cmap='RdBu',
                            norm= col.Normalize(-2, 2))
            fig.savefig(f"dH_dt_map_{bisicles_data.time:05.1f}.png")
    
    FIGS_PLOTTED = True

# Create video from the saved images

if FIGS_PLOTTED:

    print("Creating video from saved dH/dt maps...")
    image_patterns = ["./dH_dt_map_*0.png"]
    output_video_path = "./dH_dt_video.mp4"
    fps = 2  # frames per second

    create_video.create_video(output_video_path, image_patterns, fps)

