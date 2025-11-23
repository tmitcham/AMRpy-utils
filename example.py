import os
import AMRh5
import matplotlib.colors as col
from matplotlib import animation
from matplotlib import pyplot as plt
from plot_maps import plot_map, plot_diff_map, create_video

# load BISICLES h5 file
print("Loading BISICLES h5 file for thickness variable...")
bisicles_file_thickness = AMRh5.BISICLESh5(r"/mnt/c/Users/tm17544/OneDrive - University of Bristol/Projects/UKESM/UKESM1p3_initial_states/modern_AIS/hist_forcing/plot.modern.hist.000905.2d.hdf5", 'thickness')
flat_thck = bisicles_file_thickness.flatten(lev=-1)

print("Plotting thickness map...")
# plot thickness map
fig, ax = plot_map(flat_thck, 
                   kwargs_map={'title':'Ice Thickness Map',
                               'xlabel':'X (m)', 'ylabel':'Y (m)'},
                   kwargs_colorbar={'label':'Thickness (m)',
                                    'norm': col.Normalize(vmin=0, vmax=3000)}, 
                   cmap='viridis')
fig.savefig('thickness_map.png')

print("Loading another BISICLES h5 file for comparison...")
# load another BISICLES h5 file for comparison
bisicles_file_thickness_2 = AMRh5.BISICLESh5(r"/mnt/c/Users/tm17544/OneDrive - University of Bristol/Projects/UKESM/UKESM1p3_initial_states/modern_AIS/hist_forcing/plot.modern.hist.008275.2d.hdf5", 'thickness')
flat_thck_2 = bisicles_file_thickness_2.flatten(lev=-1)

print("Plotting thickness difference map...")
# plot difference map
fig_diff, ax_diff = plot_diff_map(flat_thck, flat_thck_2,
                                    kwargs_map={'title':'Thickness Difference Map',
                                               'xlabel':'X (m)', 'ylabel':'Y (m)'},
                                    kwargs_colorbar={'label':'Thickness Difference (m)',
                                                    'norm': col.Normalize(vmin=-500, vmax=500)}, 
                                    cmap='bwr')
fig_diff.savefig('thickness_difference_map.png')

# Plot every file in a directory and then create a video from the figure objects

print("Plotting thickness evolution from all files in directory...")
data_for_anim = []
data_directory = r"/mnt/c/Users/tm17544/OneDrive - University of Bristol/Projects/UKESM/UKESM1p3_initial_states/modern_AIS/hist_forcing/"
for filename in sorted(os.listdir(data_directory)):
    if filename.endswith('.hdf5') and filename.startswith('plot.'):
        print(f"Processing file: {filename}")
        bisicles_file = AMRh5.BISICLESh5(os.path.join(data_directory, filename), 'thickness').flatten(lev=-1)
        data_for_anim.append(bisicles_file)

x_coords = data_for_anim[0].x
y_coords = data_for_anim[0].y
data_frames = [obj.data for obj in data_for_anim]

def update_frame(frame_index, cax, frames, ax_obj, **kwargs_map):
        
    new_data = frames[frame_index]

    cax.set_array(new_data[:-1, :-1].ravel())

    year = data_for_anim[frame_index].time
    ax_obj.set_title(f'Thickness Map in year: {year}', **kwargs_map)
        
    return cax, ax_obj

def create_animation(data_frames, x_coords, y_coords, output_filename, fps=10, kwargs_map={}, kwargs_colorbar={}, **kwargs_pcolormesh):
    """
    Creates a video from a list of 2D data arrays using FuncAnimation.
    """
    if not data_frames:
        print("No data frames provided.")
        return

    # 1. INITIAL SETUP (Equivalent to original plot_map)
    fig, ax = plt.subplots(figsize=(8, 6))

    # Get the data for the first frame to initialize the plot
    initial_data = data_frames[0]
    
    # Initialize pcolormesh object
    # The 'animated=True' flag helps optimization
    cax_object = ax.pcolormesh(x_coords, y_coords, initial_data, shading='auto', animated=True, **kwargs_pcolormesh)
    
    # Add colorbar (only needs to be done once)
    fig.colorbar(cax_object, ax=ax, **kwargs_colorbar)
    
    # Apply static map customizations
    if 'xlabel' in kwargs_map: ax.set_xlabel(kwargs_map['xlabel'])
    if 'ylabel' in kwargs_map: ax.set_ylabel(kwargs_map['ylabel'])
    if 'xlim' in kwargs_map: ax.set_xlim(kwargs_map['xlim'])
    if 'ylim' in kwargs_map: ax.set_ylim(kwargs_map['ylim'])
    
    # 2. CREATE THE ANIMATION
    ani = animation.FuncAnimation(
        fig, 
        # The function to call for each frame
        update_frame,
        # The arguments passed to update_frame (excluding the frame index)
        fargs=(cax_object, data_frames, ax, kwargs_map),
        # The iterable defining the frames (indices 0 to N-1)
        frames=len(data_frames),
        interval=1000/fps, # Time between frames in ms
        blit=False, # Blitting is complex with pcolormesh, often safer to leave False
        repeat=False
    )
    
    # 3. SAVE THE ANIMATION
    print(f"Saving video to {output_filename}...")
    ani.save(output_filename, writer='ffmpeg', fps=fps)
    plt.close(fig)
    print("Video saved.")


create_video(figures, 'thickness_evolution.mp4', fps=5)