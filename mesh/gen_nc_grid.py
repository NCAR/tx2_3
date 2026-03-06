import numpy as np
import netCDF4
import numpy as np
from datetime import datetime

# MOM6 supergrid:
nc_sgrd = netCDF4.Dataset('../supergrid/ORCA_gridgen/ocean_hgrid_250930.nc')
x = nc_sgrd.variables['x'][:]
y = nc_sgrd.variables['y'][:]
dx = nc_sgrd.variables['dx'][:]
dy = nc_sgrd.variables['dy'][:]
area = nc_sgrd.variables['area'][:]
angle_dx = nc_sgrd.variables['angle_dx'][:]
nc_sgrd.close()

# Topography and mask
nc_topo = netCDF4.Dataset('../topography/topo.sub150.tx2_3v3.SRTM15_V2.4.edit2.SmL1.0_C1.0.nc')
depth = nc_topo.variables['D_edit2'][:]
tmask = nc_topo.variables['mask'][:]
nc_topo.close()

# T point locations
tlon = x[1::2,1::2]
tlat = y[1::2,1::2]

# U point locations
ulon = x[1::2,::2]
ulat = y[1::2,::2]

# V point locations
vlon = x[::2,1::2]
vlat = y[::2,1::2]

# Corner point locations
qlon = x[::2,::2]
qlat = y[::2,::2]

# T cell area (sum of four supergrid cells)
tarea = area[::2,::2] + area[1::2,1::2] + area[::2,1::2] + area[::2,1::2]

# t-point angle
angle = angle_dx[1::2,1::2]

# x-distance between u-points, centered at t
dxt = dx[1::2,::2] + dx[1::2,1::2]

# y-distance between v-points, centered at t
dyt = dy[::2,1::2] + dy[1::2,1::2]

# x-distance between  q-points, centered at v
dxCv = dx[2::2,::2] + dx[2::2,1::2]

# y-distance between  q-points, centered at u
dyCu = dy[::2,2::2] + dy[1::2,2::2]

# x-distance between t-points, centered at u
dxCu = dx[1::2,1::2] + np.roll(dx[1::2,1::2], -1, axis=-1)

# y-distance between t-points, centered at v
dyCv = dy[1::2,1::2] + np.roll(dy[1::2,1::2], -1, axis=0)

def write_nc_var(var, name, dimensions, long_name=None, units=None, coordinates=None):
    nc.createVariable(name, 'f8', dimensions)
    if long_name is not None:
        nc.variables[name].long_name = long_name
    if units is not None:
        nc.variables[name].units = units
    if coordinates is not None:
        nc.variables[name].coordinates = coordinates
    nc.variables[name][:] = var
    print ('... wrote variable {}'.format(name))

filename = 'tx2_3v3_grid.nc'
nc = netCDF4.Dataset(filename, 'w', format='NETCDF3_64BIT')
nc.Description = 'CESM MOM6 2/3 degree grid'
nc.Author = 'Frank, Fred, Gustavo (gmarques@ucar.edu)'
nc.Created = datetime.now().isoformat()
nc.type = 'Glogal 2/3 degree grid file'

# grid aspect ratio
ar = dyt / dxt

# grid effective grid spacing
# A = 4*pi*r^2 , area of sphere of radius r
# dA = (r*cos(theta)*dlambda)*(r*dtheta), differential area on sphere
#    = r^2*domega
# domega = dA/r^2, differential solid angle  (steradians, sr)
# 1 sr = (180./pi)^2 square degrees
costheta = np.cos(tlat*np.pi/180.)
rearth = 637122000 # Earth radius in centimeter
domega = tarea / rearth**2
egs  = np.sqrt(domega * (180./np.pi)**2)

# create netcdf dimension
M, L = qlon.shape
nc.createDimension('nyp', M)
nc.createDimension('nxp', L)
nc.createDimension('ny', M-1)
nc.createDimension('nx', L-1)


# write variable
write_nc_var(tlon, 'tlon', ('ny', 'nx'), 'array of t-grid longitudes', 'degrees_east')
write_nc_var(tlat, 'tlat', ('ny', 'nx'), 'array of t-grid latitudes', 'degrees_north')
write_nc_var(ulon, 'ulon', ('ny', 'nxp'), 'array of u-grid longitudes', 'degrees_east')
write_nc_var(ulat, 'ulat', ('ny', 'nxp'), 'array of u-grid latitudes', 'degrees_north')
write_nc_var(vlon, 'vlon', ('nyp', 'nx'), 'array of v-grid longitudes', 'degrees_east')
write_nc_var(vlat, 'vlat', ('nyp', 'nx'), 'array of v-grid latitudes', 'degrees_north')
write_nc_var(qlon, 'qlon', ('nyp', 'nxp'), 'array of q-grid longitudes', 'degrees_east')
write_nc_var(qlat, 'qlat', ('nyp', 'nxp'), 'array of q-grid latitudes', 'degrees_north')

write_nc_var(dxt, 'dxt', ('ny', 'nx'), 'x-distance between u-points, centered at t', 'meters')
write_nc_var(dyt, 'dyt', ('ny', 'nx'), 'y-distance between v-points, centered at t', 'meters')
write_nc_var(dxCv, 'dxCv', ('ny', 'nx'), 'x-distance between  q-points, centered at v', 'meters')
write_nc_var(dyCu, 'dyCu', ('ny', 'nx'), 'y-distance between  q-points, centered at u', 'meters')
write_nc_var(dxCu, 'dxCu', ('ny', 'nx'), 'x-distance between  t-points, centered at u', 'meters')
write_nc_var(dyCv, 'dyCv', ('ny', 'nx'), 'y-distance between  t-points, centered at v', 'meters')
write_nc_var(tarea, 'tarea', ('ny', 'nx'), 'area of t-cells', 'meters^2')
write_nc_var(tmask, 'tmask', ('ny', 'nx'), 'ocean fraction at t-cell centers', 'none')
write_nc_var(angle, 'angle', ('ny', 'nx'), 'angle grid makes with latitude line', 'degrees')
write_nc_var(depth, 'depth', ('ny', 'nx'), 'depth at h points', 'meters')
write_nc_var(ar, 'ar', ('ny', 'nx'), 'grid aspect ratio (dyt/dxt)', 'none')
write_nc_var(egs, 'egs', ('ny', 'nx'), 'grid effective grid spacing', 'degrees')

# close files
nc.close()





