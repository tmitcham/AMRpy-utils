import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator

class AMRh5:

    def __init__(self, filename, variable_name):

        self.file = None
        self.filename = filename
        self.variable_name = variable_name
        self.time = None
        self.num_levels = None
        self.offsets = None
        self.num_boxes = None
        self.dx = None
        self.x = None
        self.y = None
        self.data = None

        self.open()

        n_components = self.file.attrs['num_components']

        self.find_var(n_components)

        target_component = self.find_var(n_components)

        self.read_var(target_component)

        self.close()

    def open(self):

        print(f"Opening file: {self.filename}")
        self.file = h5py.File(self.filename, 'r')
        if self.file is None:
            raise RuntimeError(f"Could not open file: {self.filename}")
        
    def close(self):

        if self.file is not None:
            self.file.close()
            self.file = None

    def find_var(self, n_components):

        count = 0

        while count < n_components and self.file.attrs['component_'+str(count)].decode('utf-8') != self.variable_name:
            count += 1

        if count == n_components:
            count = np.nan
            raise KeyError(f"Variable '{self.variable_name}' not found in the HDF5 file.")
        
        return count

    def read_var(self, target_component):

        if self.file is None:
            raise RuntimeError("HDF5 file is not opened.")
        
        print(f"Reading variable '{self.variable_name}' at component index: {target_component}")

        n_level = self.file.attrs['num_levels']

        # create stubs for various quantities that will be read these will be
        # constructed using .append
        levelsdx = [None] * n_level
        levelsboxes = [None] * n_level
        levelsoff = [None] * n_level
        levelsx = [None] * n_level
        levelsy = [None] * n_level
        levelsdata = [None] * n_level

        # crrate a time variable in case it isnt picked up from level info
        time = np.nan

        # visit each level in turn coarse to fine
        for level in range(n_level):

            print(f"Reading data on level {level}: Level {level+1}/{n_level}")

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

            # print(dx,offset)

            # read model time if it is there (years)
            if 'time' in self.file.attrs:
                time = h5level.attrs['time']

            # read data (h5data) for this level, which is found in boxes so need to read
            # pointer into data for each box (h5offs)

            # print(h5level.keys())
            # ['Processors', 'boxes', 'data:datatype=0', 'data:offsets=0', 'data_attributes']
            h5box = np.array(h5level.get('boxes'))
            h5data = np.array(h5level.get('data:datatype=0'))
            h5offs = np.array(h5level.get('data:offsets=0'))

            n_boxes = len(h5box)

            # create stubs for quantities (x,y,data) held in boxes on this level
            boxesx = [None] * n_boxes
            boxesy = [None] * n_boxes
            boxesdata = [None] * n_boxes

            # visit each box in this level
            for box in range(n_boxes):

                # find x, y for this box from bounding box information remembering
                # level offsets in x,y and add border of ghost cells
                # NOT clear why +2 but provides correct numbers when compared to h5offs
                y = np.arange(h5box['lo_i'][box]-1,h5box['hi_i'][box]+2) * dx + offset
                x = np.arange(h5box['lo_j'][box]-1,h5box['hi_j'][box]+2) * dx + offset

                # add x, y data to evolving level
                boxesx[box] = x
                boxesy[box] = y

                # box size
                nx = len(x)
                ny = len(y)

                # find locations of box data in h5
                en = h5offs[box]

                # if level < 2:
                #   print(level, box, h5offs[box], h5box['lo_i'][box], h5box['hi_i'][box], h5box['lo_j'][box], h5box['hi_j'][box], nx, ny)

                st = en + target_component * nx * ny
                en = st + nx * ny

                boxesdata[box] = h5data[st:en].reshape((nx,ny))

            # once completed all boxes in a level add them to the overall structure
            levelsx[level] = boxesx
            levelsy[level] = boxesy
            levelsdata[level] = boxesdata

            # also add useful level info such as grid spacing, grid offset and number boxes
            levelsdx[level] = dx
            levelsoff[level] = offset
            levelsboxes[level] = n_boxes

        # pack all of the information into the class
        self.time = time
        self.num_levels = n_level # number of levels
        self.offset = levelsoff # grid offsets for each level
        self.num_boxes = levelsboxes # list of number of boxes in each level
        self.dx = levelsdx # grid spacing for each level
        self.x = levelsx # x, y for each level in global coordinate system
        self.y = levelsy # as list of arrays [level][box]
        self.data = levelsdata # data for each level as list of arrays [level][box]

    def flatten(self,lev=-1,xmin=np.nan,xmax=np.nan,ymin=np.nan,ymax=np.nan):

        print(f"Flattening data to level {lev}")

        return flatAMRh5(self,lev,xmin,xmax,ymin,ymax)
    
class flatAMRh5:

    def __init__(self, AMRh5Obj, target_level=-1, xmin=np.nan, xmax=np.nan, ymin=np.nan, ymax=np.nan):

        self.fname = AMRh5Obj.filename
        self.time = AMRh5Obj.time
        self.variable_name = AMRh5Obj.variable_name
        self.dx = AMRh5Obj.dx[target_level]
        self.x = None
        self.y = None
        self.data = None

        # if bounding box not specfied then
        # find based on x,y bounds of the coarsest grid
        # the reason for doing + dx0 on the min values is due to a single ghost cell layer

        dx0 = AMRh5Obj.dx[0]/2.0
        xmin = np.min(AMRh5Obj.x[0][:]) + dx0
        xmax = np.max(AMRh5Obj.x[0][:]) - dx0
        ymin = np.min(AMRh5Obj.y[0][:]) + dx0
        ymax = np.max(AMRh5Obj.y[0][:]) - dx0

        print(xmin,xmax,ymin,ymax)

        # create x,y coordinates for base grid at resolution specfied in call (default is finest)
        # xx = np.arange(xmin,xmax+self.dx[lev],self.dx[lev]) + self.offset[lev]
        # yy = np.arange(ymin,ymax+self.dx[lev],self.dx[lev]) + self.offset[lev]

        # stopping at xmax/ymax ensures that the final value is always less than the max value due to less than in arange

        xx = np.arange(xmin+0.5*AMRh5Obj.dx[target_level],xmax,AMRh5Obj.dx[target_level])
        yy = np.arange(ymin+0.5*AMRh5Obj.dx[target_level],ymax,AMRh5Obj.dx[target_level])

        # if bounds have been given find nearest poinst on base grid
        if not np.isnan(xmin):
            xx = xx[np.nonzero(xx>=xmin)]

        if not np.isnan(xmax):
            xx = xx[np.nonzero(xx<=xmax)]

        if not np.isnan(ymin):
            yy = yy[np.nonzero(yy>=ymin)]

        if not np.isnan(ymax):
            yy = yy[np.nonzero(yy<=ymax)]

        print(f"Flattened grid size: {len(xx)} x {len(yy)}")
        print(f"Bounding box for flattened data: xmin={xx[0]}, xmax={xx[-1]}, ymin={yy[0]}, ymax={yy[-1]}")
        print(f"Grid spacing for flattened data: dx={xx[1]-xx[0]}")

        # create suitably sized array to hold flattened data (base grid)
        data = np.zeros((len(xx),len(yy)))

        # visit each level (coarsest upwards) until you get to the target level, 
        # and then visit each box
        for level in range(AMRh5Obj.num_levels if target_level < 0 else target_level+1):

            print(f"Interpolating level {level} onto target grid. Level {level+1}/{AMRh5Obj.num_levels}")

            for box in range(AMRh5Obj.num_boxes[level]):

                # find the i,j corrdinates of the bounds for this box
                # on the new base grid

                imin = self.findend(xx,AMRh5Obj.x[level][box],0,1)
                imax = self.findend(xx,AMRh5Obj.x[level][box],-2,-1) + 2

                jmin = self.findend(yy,AMRh5Obj.y[level][box],0,1)
                jmax = self.findend(yy,AMRh5Obj.y[level][box],-2,-1) + 2

                # create and interpolator for this box
                # including any ghost cells
                interp = RegularGridInterpolator((AMRh5Obj.x[level][box],AMRh5Obj.y[level][box]), \
                                                AMRh5Obj.data[level][box], \
                                                method='nearest',bounds_error=False,fill_value=None)

                # create x,y coordinates for poiints on new base grid that fall within this box
                # and are to be filled
                xm,ym = np.meshgrid(yy[jmin:jmax],xx[imin:imax])

                # print(imin,imax,jmin,jmax)
                # print(np.shape(yy[jmin:jmax]),np.shape(xx[imin:imax]))
                # print(np.shape(xm),np.shape(ym),np.shape(data[imin:imax,jmin:jmax]))

                # do the nearest neighbor interpolation
                data[imin:imax,jmin:jmax] = interp((ym,xm))

        self.x = xx
        self.y = yy
        self.data = data

    def findend(self,xx,x,i1,i2):

        value = np.where((xx >= x[i1]) & (xx <= x[i2]))

        return value[0][0]
    
class BISICLESh5(AMRh5):

    def __init__(self, filename, variable_name):

        super().__init__(filename, variable_name)

        self.full_name = None
        self.units = None
        self.bisicles_h5_attrs = {}

        print(f"Getting full name and units for variable '{self.variable_name}'.")

        self.full_name, self.units = self.get_full_name_units(self.variable_name)

        self.open()

        print(f"Reading plotfile header attributes from file: {self.filename}")

        self.bisicles_h5_attrs = self.read_plotfile_attrs()

        self.close()

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
        
    
    def read_plotfile_attrs(self):

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

            if key in self.file.attrs:

                if isinstance(self.file.attrs[key], bytes):
                    attrs[key] = self.file.attrs[key].decode('utf-8')
                else:
                    attrs[key] = self.file.attrs[key].item()
        
        return attrs
    
    def flatten(self,lev=-1,xmin=np.nan,xmax=np.nan,ymin=np.nan,ymax=np.nan):

        print(f"Flattening data to level {lev}")

        return flatBISICLESh5(self,lev,xmin,xmax,ymin,ymax)


class flatBISICLESh5(flatAMRh5):

    def __init__(self, BISICLESh5Obj, target_level=-1, xmin=np.nan, xmax=np.nan, ymin=np.nan, ymax=np.nan):

        super().__init__(BISICLESh5Obj, target_level, xmin, xmax, ymin, ymax)

        self.full_name = BISICLESh5Obj.full_name
        self.units = BISICLESh5Obj.units
        self.bisicles_h5_attrs = BISICLESh5Obj.bisicles_h5_attrs