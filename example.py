import os
import AMRh5
import matplotlib.colors as col
from plot_maps import plot_map, plot_diff_map, create_video

# load BISICLES h5 file
bisicles_file_thickness = AMRh5.BISICLESh5('path_to_bisicles_file.h', 'thickness')
flat_thck = bisicles_file_thickness.flatten()

# plot thickness map
fig, ax = plot_map(flat_thck, 
                   kwargs_map={'title':'Ice Thickness Map',
                               'xlabel':'X (m)', 'ylabel':'Y (m)'},
                   kwargs_colorbar={'label':'Thickness (m)',
                                    'norm': col.Normalize(vmin=0, vmax=3000)}, 
                   cmap='viridis')
fig.savefig('thickness_map.png')

# load another BISICLES h5 file for comparison
bisicles_file_thickness_2 = AMRh5.BISICLESh5('path_to_another_bisicles_file.h', 'thickness')
flat_thck_2 = bisicles_file_thickness_2.flatten()

# plot difference map
fig_diff, ax_diff = plot_diff_map(flat_thck, flat_thck_2,
                                    kwargs_map={'title':'Thickness Difference Map',
                                               'xlabel':'X (m)', 'ylabel':'Y (m)'},
                                    kwargs_colorbar={'label':'Thickness Difference (m)',
                                                    'norm': col.Normalize(vmin=-500, vmax=500)}, 
                                    cmap='bwr')
fig_diff.savefig('thickness_difference_map.png')

# Plot every file in a directory and then create a video from the figure objects

figures = []
data_directory = 'path_to_bisicles_files_directory'
for filename in sorted(os.listdir(data_directory)):
    if filename.endswith('.hdf5'):
        bisicles_file = AMRh5.BISICLESh5(os.path.join(data_directory, filename), 'thickness')
        flat_thck = bisicles_file.flatten()
        year = flat_thck.time
        
        fig, ax = plot_map(flat_thck, 
                           kwargs_map={'title':f'Thickness Map in year: {year}',        
                                       'xlabel':'X (m)', 'ylabel':'Y (m)'},
                           kwargs_colorbar={'label':'Thickness (m)',
                                            'norm': col.Normalize(vmin=0, vmax=3000)}, 
                           cmap='viridis')
        figures.append(fig)

create_video(figures, 'thickness_evolution.mp4', fps=5)