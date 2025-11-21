# Take an AMRh5 object and flatten the data into a uniform grid

from AMRh5 import AMRh5, bisiclesh5

import numpy as np
from scipy.interpolate import RegularGridInterpolator

class AMRflatten:

    def __init__(self, amrh5, lev=-1, xmin=np.nan, xmax=np.nan, ymin=np.nan, ymax=np.nan):

        # if bounding box not specfied then
        # find based on x,y bounds of the coarsest grid
        # if bounds have been given find nearest poinst on base grid
        # create x,y coordinates for base grid at resolution specfied in call (default is finest)

        # Why is this not the calculation we use?
        # xx = np.arange(xmin,xmax+self.dx[lev],self.dx[lev]) + self.offset[lev]
        # yy = np.arange(ymin,ymax+self.dx[lev],self.dx[lev]) + self.offset[lev]

        dx0 = amrh5.dx[0]/2.0
        xmin_default = np.min(amrh5.x[0][:]) + dx0
        xmax_default = np.max(amrh5.x[0][:]) - dx0
        ymin_default = np.min(amrh5.y[0][:]) + dx0
        ymax_default = np.max(amrh5.y[0][:]) - dx0

        xx = np.arange(xmin_default+0.5*amrh5.dx[lev],xmax_default-1.5*amrh5.dx[lev],amrh5.dx[lev])
        yy = np.arange(ymin_default+0.5*amrh5.dx[lev],ymax_default-1.5*amrh5.dx[lev],amrh5.dx[lev])

        if np.isnan(xmin):
            xmin = np.min(amrh5.x[0][:]) + dx0
        else:
            xx = np.arange(xmin+0.5*amrh5.dx[lev],xmax-1.5*amrh5.dx[lev],amrh5.dx[lev])
            xmin = xx[np.nonzero(xx>xmin)[0][[0]]]

        if np.isnan(xmax):
            xmax = np.max(amrh5.x[0][:]) - dx0
        else:
            xx = np.arange(xmin+0.5*amrh5.dx[lev],xmax-1.5*amrh5.dx[lev],amrh5.dx[lev])
            xmax = xx[np.nonzero(xx<xmax)[0][[-1]]]

        if np.isnan(ymin):
            ymin = np.min(amrh5.y[0][:]) + dx0
        else:
            yy = np.arange(ymin+0.5*amrh5.dx[lev],ymax-1.5*amrh5.dx[lev],amrh5.dx[lev])
            ymin = yy[np.nonzero(yy>ymin)[0][[0]]]

        if np.isnan(ymax):
            ymax = np.max(amrh5.y[0][:]) - dx0
        else:
            yy = np.arange(ymin+0.5*amrh5.dx[lev],ymax-1.5*amrh5.dx[lev],amrh5.dx[lev])
            ymax = yy[np.nonzero(yy<ymax)[0][[-1]]]
        
        # for testing
        xx = np.arange(np.min(amrh5.x[0][:])+dx0,np.max(amrh5.x[0][:])-1.5*amrh5.dx[0],amrh5.dx[0])
        yy = np.arange(np.min(amrh5.y[0][:])+dx0,np.max(amrh5.y[0][:])-1.5*amrh5.dx[0],amrh5.dx[0])
        print(xmin,xmax,ymin,ymax)
        print(xx[0],xx[1],xx[-2],xx[-1])
        print(yy[0],yy[1],yy[-2],yy[-1])

        # create suitably sized array to hold flattened data (base grid)
        data = np.zeros((len(xx),len(yy)))

    # visit each level (coarest downwards) and box

    def findend(xx,x,i1,i2):

      value = np.where((xx >= x[i1]) & (xx <= x[i2]))

      return value[0][0]
    
    
    def flatten_data(self):
      for level in range(self.num_levels):
        for box in range(self.num_boxes[level]):

            # print(level,box)

            # find the i,j corrdinates of the bounds for this box
            # on the new base grid

            imin = findend(xx,self.x[level][box],0,1)
            imax = findend(xx,self.x[level][box],-2,-1) + 2

            jmin = findend(yy,self.y[level][box],0,1)
            jmax = findend(yy,self.y[level][box],-2,-1) + 2

            # create and interpolator for this box
            # including any ghost cells
            interp = RegularGridInterpolator((self.x[level][box],self.y[level][box]), \
                                            self.data[level][box], \
                                            method='nearest',bounds_error=False,fill_value=None)

            # create x,y coordinates for poiints on new base grid that fall within this box
            # and are to be filled
            xm,ym = np.meshgrid(yy[jmin:jmax],xx[imin:imax])

            # print(imin,imax,jmin,jmax)
            # print(np.shape(yy[jmin:jmax]),np.shape(xx[imin:imax]))
            # print(np.shape(xm),np.shape(ym),np.shape(data[imin:imax,jmin:jmax]))

            # do the nearest neighbor interpolation
            data[imin:imax,jmin:jmax] = interp((ym,xm))


    
    

    

    print(xmin,xmax,ymin,ymax)
    print(xx[0],xx[1],xx[-2],xx[-1])
    print(yy[0],yy[1],yy[-2],yy[-1])

    # create suitably sized array to hold flattened data (base grid)
    data = np.zeros((len(xx),len(yy)))

    # visit each level (coarest downwards) and box
    for level in range(self.num_levels):
      for box in range(self.num_boxes[level]):

        # print(level,box)

        # find the i,j corrdinates of the bounds for this box
        # on the new base grid

        imin = findend(xx,self.x[level][box],0,1)
        imax = findend(xx,self.x[level][box],-2,-1) + 2

        jmin = findend(yy,self.y[level][box],0,1)
        jmax = findend(yy,self.y[level][box],-2,-1) + 2

        # create and interpolator for this box
        # including any ghost cells
        interp = RegularGridInterpolator((self.x[level][box],self.y[level][box]), \
                                         self.data[level][box], \
                                         method='nearest',bounds_error=False,fill_value=None)

        # create x,y coordinates for poiints on new base grid that fall within this box
        # and are to be filled
        xm,ym = np.meshgrid(yy[jmin:jmax],xx[imin:imax])

        # print(imin,imax,jmin,jmax)
        # print(np.shape(yy[jmin:jmax]),np.shape(xx[imin:imax]))
        # print(np.shape(xm),np.shape(ym),np.shape(data[imin:imax,jmin:jmax]))

        # do the nearest neighbor interpolation
        data[imin:imax,jmin:jmax] = interp((ym,xm))

        # create a class to hold flattened version and return it

    return flat(self.fname,self.time,self.variable_name,self.full_name,self.units,self.dx[lev],xx,yy,data)