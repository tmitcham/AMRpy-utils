import h5py
import numpy as np

class AMRh5:

    def __init__(self, filename, variable_name, read_data=True, write_data=False):
        
        self.file = None
        self.file_name = filename
        self.time = None
        self.variable_name = variable_name
        self.num_levels = None # number of levels
        self.offset = None # grid offsets for each level
        self.num_boxes = None # list of number of boxes in each level
        self.dx = None # grid spacing for each level
        self.x = None # x, y for each level in global coordinate system
        self.y = None # as list of arrays [level][box]
        self.data = None # data for each level as list of arrays [level][box]

        self.open(filename)
        var_comp_num = self.find_var(variable_name)

        if read_data:
            self.read_var(var_comp_num)

        if write_data:
            pass # to be implemented

    def open(self, filename):
        self.file = h5py.File(filename, 'r')
        self.file_name = filename

    def close(self):
        if self.file is not None:
            self.file.close()
            self.file = None
    
    def find_var(self, variable_name):

        if self.file is None:   
            raise RuntimeError("HDF5 file is not opened.")
        
        n_components = self.file.attrs['num_components']
        
        count = 0

        while count < n_components and self.file.attrs['component_'+str(count)].decode('utf-8') != variable_name:
            count += 1

        if count == n_components:
            raise KeyError(f"Variable '{variable_name}' not found in the HDF5 file.")
        
        self.variable_name = variable_name

        return count
    
    def read_var(self, var_comp_num):

        if self.file is None:
            raise RuntimeError("HDF5 file is not opened.")
        
        n_levels = self.file.attrs['num_levels']

        # create stubs for various quantities that will be read these will be
        # constructed using .append
        levelsdx = []
        levelsboxes = []
        levelsoff = []
        levelsx = []
        levelsy = []
        levelsdata = []

        # crrate a time variable in case it isnt picked up from level info
        time = np.nan

        # visit each level in turn coarse to fine
        for level in range(n_levels):
        
            h5level = self.file['/level_'+str(level)+'/']

            # h5level.attrs.keys()
            # ['dt', 'dx', 'prob_domain', 'ref_ratio', 'time']

            # read this levels grid spacing and calculate its offset
            dx = h5level.attrs['dx']

            # this is the offset in x,y required to align grids at different levels
            if level == 0:
                offset = dx / 2.0
            else:
                offset = offset - dx / 2.0

            # read model time if it is there (years)
            if h5level.attrs.__contains__('time'):
                time = h5level.attrs['time']

            # read data (h5data) for this level, which is found in boxes so need to read
            # pointer into data for each box (h5offs)

            h5box = h5level.get('boxes')
            h5box = np.array(h5box)
            h5data = h5level.get('data:datatype=0')
            h5data = np.array(h5data)
            h5offs = h5level.get('data:offsets=0')
            h5offs = np.array(h5offs)

            n_boxes = len(h5box)

            # create stubs for quantities (x,y,data) held in boxes on this level
            boxesx = []
            boxesy = []
            boxesdata = []

            # visit each box in this level
            for box in range(n_boxes):

                # find x, y for this box from bounding box information remembering
                # level offsets in x,y and add border of ghost cells
                # NOT clear why +2 but provides correct numbers when compared to h5offs (think it's because of python less than indexing)
                y = np.arange(h5box['lo_i'][box]-1,h5box['hi_i'][box]+2) * dx + offset
                x = np.arange(h5box['lo_j'][box]-1,h5box['hi_j'][box]+2) * dx + offset

                # add x, y data to evolving level
                boxesx.append(x)
                boxesy.append(y)

                # box size
                nx = len(x)
                ny = len(y)

                # find locations of box data in h5
                en = h5offs[box]

                # if level < 2:
                # print(level, box, h5offs[box], h5box['lo_i'][box], h5box['hi_i'][box], h5box['lo_j'][box], h5box['hi_j'][box], nx, ny)

                st = en + var_comp_num * nx * ny
                en = st + nx * ny

                boxesdata.append(h5data[st:en].reshape((nx,ny)))

            # once completed all boxes in a level add them to the overall structure
            levelsx.append(boxesx)
            levelsy.append(boxesy)
            levelsdata.append(boxesdata)

            # also add useful level info such as grid spacing, grid offset and number boxes
            levelsdx.append(dx)
            levelsoff.append(offset)
            levelsboxes.append(n_boxes)

        # once all levels are read assign to class variables
        self.time = time
        self.num_levels = n_levels
        self.dx = levelsdx
        self.offset = levelsoff
        self.num_boxes = levelsboxes
        self.x = levelsx
        self.y = levelsy
        self.data = levelsdata

class bisiclesh5(AMRh5):

    def __init__(self, filename, variable_name):

        super().__init__(filename, variable_name)

        self.full_name = None
        self.units = None
        self.bisicles_h5_attrs = {}

        self.full_name, self.units = self.get_full_name_units(variable_name)

        self.open(filename)
        self.bisicles_h5_attrs = self.read_plotfile_attrs(filename)

    def get_full_name_units(self, vname):
    
        variable_unit_table = {
            'thickness':['Ice thickness','m'],
            'xVel':['Horizontal velocity x-component','m/yr'],
            'yVel':['Horizontal velocity y-component','m/yr'],
            'vVel':['Vertical velocity','m/yr'],
            'Z_surface':['Upper surface elevation','m'],
            'Z_bottom':['Lower surface elevation','m'],
            'Z_base':['Bedrock elevation','m'],
            'basal_friction':['Basal friction','Pa'],
            'C0':['-','-'],
            'xRhs':['Gravitational driving stress in x','Pa'],
            'yRhs':['Gravitational driving stress in y','Pa'],
            'dThickness/dt':['Rate of thickness change','m/yr'],
            'xfVel':['Horizontal velocity x-component at upper surface','m/yr'],
            'yfVel':['Horizontal velocity y-component at upper surface','m/yr'],
            'zfVel':['Vertical velocity at upper surface','m/yr'],
            'xbVel':['Horizontal velocity x-component at lower surface','m/yr'],
            'ybVel':['Horizontal velocity y-component at lower surface','m/yr'],
            'zbVel':['Vertical velocity at lower surface','m/yr'],
            'dragCoef':['Basal friction coefficient','(units depend on law)'],
            'viscosityCoef':['viscosityCoef','-'],
            'xxViscousTensor':['xxViscousTensor','Pa'],
            'yxViscousTensor':['yxViscousTensor','Pa'],
            'xyViscousTensor':['xyViscousTensor','Pa'],
            'yyViscousTensor':['yyViscousTensor','Pa'],
            'activeBasalThicknessSource':['Mass balance at lower surface (active)','m/yr'],
            'activeSurfaceThicknessSource':['Mass balance at upper surface (active)','m/yr'],
            'divergenceThicknessFlux':['Divergence of horizontal flux','m/yr'],
            'basalThicknessSource':['Mass balance at lower surface','m/yr'],
            'surfaceThicknessSource':['Mass balance at upper surface','m/yr'],
            'calvingFlux':['Calving flux','-'],
            'mask':['mask','-'] 
        }

        if vname in variable_unit_table:
            return variable_unit_table[vname][0], variable_unit_table[vname][1]
        
    
    def read_plotfile_attrs(self, filename):

        if self.file is None:
            raise RuntimeError("HDF5 file is not opened.")
        
        keys = ['bisicles_patch_number', 
                'bisicles_version_major', 
                'bisicles_version_minor', 
                'chombo_patch_number', 
                'chombo_version_major', 
                'chombo_version_minor', 
                'crs_origin_x', 
                'crs_origin_y', 
                'current_step', 
                'density_of_ice', 
                'dt', 
                'filetype', 
                'finest_level', 
                'git_hash', 
                'git_remote', 
                'max_level', 
                'num_components', 
                'num_levels', 
                'seconds_per_unit_time', 
                'svn_repository', 
                'svn_url', 
                'svn_version', 
                'time'
                ]
        
        attrs = {}

        for key in keys:

            if self.file.attrs.__contains__(key):

                if isinstance(self.file.attrs[key], bytes):
                    attrs[key] = self.file.attrs[key].decode('utf-8')
                else:
                    attrs[key] = self.file.attrs[key].item()
        
        self.bisicles_h5_attrs = attrs

            