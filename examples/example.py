"""
Example script demonstrating the use of AMRpy-utils to load BISICLES h5 files,
flatten data, and plot thickness maps and differences with grounding lines
"""

import os
import matplotlib.colors as col
import amrpy_utils
from plot_maps import plot_map, plot_diff_map, plot_gls

# load BISICLES h5 file
print("Loading BISICLES h5 file for thickness variable...")
#thickness = AMRh5.BISICLESh5("/mnt/c/Users/tm17544/OneDrive - University of Bristol/Projects/UKESM/UKESM1p3_initial_states/modern_AIS/hist_forcing/plot.modern.hist.000905.2d.hdf5", 'thickness')
thickness = AMRh5.BISICLESh5("C:\\Users\\tm17544\\OneDrive - University of Bristol\\Projects\\UKESM\\UKESM1p3_initial_states\\modern_AIS\\hist_forcing\\plot.modern.hist.000905.2d.hdf5", 'thickness')
flat_thickness = thickness.flatten(lev=1)

#thickness_2 = AMRh5.BISICLESh5("/mnt/c/Users/tm17544/OneDrive - University of Bristol/Projects/UKESM/UKESM1p3_initial_states/modern_AIS/hist_forcing/plot.modern.hist.008275.2d.hdf5", 'thickness')
thickness_2 = AMRh5.BISICLESh5("C:\\Users\\tm17544\\OneDrive - University of Bristol\\Projects\\UKESM\\UKESM1p3_initial_states\\modern_AIS\\hist_forcing\\plot.modern.hist.008275.2d.hdf5", 'thickness')
flat_thickness_2 = thickness_2.flatten(lev=1)

#bedrock = AMRh5.BISICLESh5("/mnt/c/Users/tm17544/OneDrive - University of Bristol/Projects/UKESM/UKESM1p3_initial_states/modern_AIS/hist_forcing/plot.modern.hist.000905.2d.hdf5", 'Z_base')
bedrock = AMRh5.BISICLESh5("C:\\Users\\tm17544\\OneDrive - University of Bristol\\Projects\\UKESM\\UKESM1p3_initial_states\\modern_AIS\\hist_forcing\\plot.modern.hist.000905.2d.hdf5", 'Z_base')
flat_bedrock = bedrock.flatten(lev=1)

print("Plotting thickness map...")
# plot thickness map
fig, ax = plot_map(flat_thickness,
                   kwargs_map={'title':f"Thickness in year: {flat_thickness.bisicles_h5_attrs['time']}",
                               'xlabel':'X (m)', 'ylabel':'Y (m)'},
                   kwargs_colorbar={'label':'Thickness (m)'},
                   kwargs_pcolormesh={'cmap':'Blues',
                                      'norm': col.Normalize(0, 1000)})

contour = plot_gls(fig, ax, [flat_bedrock, flat_thickness, flat_thickness_2],
                   kwargs_contour={'colors':['black', 'magenta'],
                                      'linewidths':[1., 1.]})
fig.savefig('thickness_gl_test.png')
"""

# load BISICLES h5 file
print("Loading BISICLES h5 file for thickness variable...")
bisicles_file_1 = AMRh5.BISICLESh5("/mnt/c/Users/tm17544/OneDrive - University of Bristol/Projects/UKESM/UKESM1p3_initial_states/pi_AIS/pi_forcing/plot.pi.pi.003256.2d.hdf5", 'thickness')
flat_thickness = bisicles_file_1.flatten(lev=0)

bisicles_file_2 = AMRh5.BISICLESh5("/mnt/c/Users/tm17544/OneDrive - University of Bristol/Projects/UKESM/UKESM1p3_initial_states/pi_AIS/pi_forcing/plot.pi.pi.008407.2d.hdf5", 'thickness')
flat_thickness_2 = bisicles_file_2.flatten(lev=0)

print("Plotting thickness map...")
# plot thickness map
fig, ax = plot_diff_map(flat_thickness, flat_thickness_2,
                   kwargs_map={'title':f"Cumulative change in thickness: {flat_thickness_2.time - flat_thickness.time} years",
                               'xlabel':'X (m)', 'ylabel':'Y (m)'},
                   kwargs_colorbar={'label':'Thickness change (m)'},
                   cmap='RdBu',
                   norm= col.Normalize(-100, 100))
fig.savefig('thickness_change_map_pi_pi.png')
"""
# Plot every file in a directory and then create a video from the figure objects

if False:
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
