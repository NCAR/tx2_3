import numpy as np
import xarray as xr

import matplotlib.pyplot as plt

import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point

def mom6_point(ds,lon,lat):
    jplt = (np.abs(ds.geolat[:,0].data[:])-lat).argmin()
    for n in range(3) :
        iplt = (np.abs((ds.geolon[jplt,:].data[:]) - lon)).argmin()
        jplt = (np.abs(ds.geolat[:,iplt].data[:] - lat)).argmin()
#        print(iplt,jplt,ds.geolat[jplt,iplt].values,ds.geolon[jplt,iplt].values)
    return iplt,jplt

def mom6_latlon2ij(lon2D,lat2D,lon,lat,max_iter=10):
    ilast=-1
    jlast=-1
       
    # Find the equator
    jpt = np.abs(lat2D[:,0]).argmin().values
    
    # Handle the wierd MOM6 grids that go past -180
    lo = 0.
    if ( lon < lon2D[jpt,:].min() ) :
        lo = 360.
    if ( lon > lon2D[jpt,:].max() ) :
        lo = -360.
    # Search of the point along grid directions
    ipt = np.abs(lon2D[jpt,:] - (lon+lo)).argmin().values
    n=0
    while ( (ilast != ipt ) and (jlast != jpt) and (n < max_iter) ):

        lo = 0.
        if ( lon < lon2D[jpt,:].min() ) :
            lo = 360.
        if ( lon > lon2D[jpt,:].max() ) :
            lo = -360.
   
        ilast = ipt
        jlast = jpt
        jpt = np.abs(lat2D[:,ipt]-lat).argmin().values
        ipt = np.abs(lon2D[jpt,:] - (lon+lo)).argmin().values
        n = n+1
#        print('n=',n,', j,i=',jpt,ipt,
#              ' lon=',lon2D[jpt,ipt].values,
#              ' lat=',lat2D[jpt,ipt].values,'lo=',lo)
    return ipt,jpt
    
def map_mask(df,vmask,lonbeg,lonend,latbeg,latend,
             label_skip=10,line_skip=1,           
             vlon='geolon',vlone='geolonb',
             vlat='geolat',vlate='geolatb',
             ):
    mask = df[vmask]
    lonc = df[vlon]
    latc = df[vlat]
    lone = df[vlone]
    late = df[vlate]
    
    nlath,nlonh=mask.shape
#    print('nlonh=',nlonh,'nlath=',nlath)
    
#   Find the i,j of the corners of the domain
    ivals=np.zeros(4,dtype='int')
    jvals=np.zeros(4,dtype='int')
    ivals[0],jvals[0] = mom6_latlon2ij(lonc,latc,lonbeg,latbeg)
    ivals[1],jvals[1] = mom6_latlon2ij(lonc,latc,lonend,latbeg)
    ivals[2],jvals[2] = mom6_latlon2ij(lonc,latc,lonend,latend)
    ivals[3],jvals[3] = mom6_latlon2ij(lonc,latc,lonbeg,latend)

#   Find the extrema of the corners
    imin = max([ivals.min()-1,0])
    jmin = max([jvals.min()-1,0])
    imax = min([ivals.max()+1,nlonh])
    jmax = min([jvals.max()+1,nlath])

#   The search above has trouble near the tripole boundary
#   punt and use the whole grid
    if ( latend > 66. ) :
        imin=0
        imax = nlonh
        jmax = nlath
    
#    print('imin=',imin,' imax=',imax,' jmin=',jmin,' jmax=',jmax)

    if ( (imin >= imax) or (jmin >= jmax) ):
        print('null region imin,imax=',imin,imax,' jmin,jmax=',jmin,jmax)
        return
    
    fig,ax=plt.subplots(figsize=(10,7.5),constrained_layout=True,
                        subplot_kw={'projection': ccrs.PlateCarree()})
    plt.set_cmap('bwr_r')
    
    lo=0
    if ( lonbeg < -180. ):
        lo=360.

    pc=ax.pcolormesh(lonc[jmin:jmax,imin:imax]+lo,
                     latc[jmin:jmax,imin:imax],
                     mask[jmin:jmax,imin:imax],vmin=0,vmax=1)
    cc=ax.contour(lonc[jmin:jmax,imin:imax]+lo,
                  latc[jmin:jmax,imin:imax],
                  mask[jmin:jmax,imin:imax],
               levels=[0.5],colors='black',linewidths=1.5,linestyles='dashed')
#    ax.set_title('Unmodified Topography')
    ax.set_xlim(lonbeg+lo,lonend+lo)
    ax.set_ylim(latbeg,latend)
    
    # Longitude labels
    ax.set_xticks(np.arange(lonbeg+lo,lonend+lo), crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)

    # Latitude labels
    ax.set_yticks(np.arange(latbeg,latend), crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    ax.yaxis.set_major_formatter(lat_formatter)

    ax.coastlines()
    ax.add_feature(cfeature.COASTLINE, color='darkgrey')
    
    for i in range(imin,imax,line_skip):
        ax.plot(lone[:,i]+lo,late[:,i],
                color='grey',linewidth=0.5)
    for j in range(jmin,jmax,line_skip):
        ax.plot(lone[j,:]+lo,late[j,:],
                color='grey',linewidth=0.5)

    cbar = fig.colorbar(pc,shrink=0.85)
    for j in range(jmin,jmax,label_skip):
        for i in range(imin,imax,label_skip):
            if ( latbeg < latc[j,i] and
                 latend > latc[j,i] and
                 lonbeg+lo < lonc[j,i]+lo and
                 lonend+lo > lonc[j,i]+lo ) :
                       
                lab = f"({j},{i})"
                ax.text(lonc[j,i]+lo,latc[j,i],
                        lab,ha='center',va='center',fontsize=8)
                
def map_maskij(df,vmask,imin,imax,jmin,jmax,
               label_skip=10,line_skip=1,draw_coastline=False):
    mask = df[vmask]
    nx = df[vmask].shape[1]
    ny = df[vmask].shape[0]
    
    if ( (imin >= imax) or (jmin >= jmax) ):
        print('null region imin,imax=',imin,imax,' jmin,jmax=',jmin,jmax)
        return

    x = np.arange(imin,imax)
    y = np.arange(jmin,jmax)
    
    if  imin < 0 :
        print('imin= ',imin,"  < 0 . map_maskij can't span cyclic boundary")
        return
    elif imax > nx  :
        print('imax = ',imax,'  > ',nx," . map_maskij can't span cyclic boundary")
        return
    if  jmin < 0 :
        print('jmin= ',jmin,"  < 0 . map_maskij can't span boundary")
        return
    elif jmax > ny  :
        print('jmax = ',jmax,'  > ',ny," . map_maskij can't span boundary")
        return
              
    c = mask[jmin:jmax,imin:imax]

    fig,ax=plt.subplots(figsize=(10,7.5),constrained_layout=True)
    plt.set_cmap('bwr_r')
    
    pc=ax.pcolormesh(x,y,c,vmin=0,vmax=1)
#    ax.set_title('Unmodified Topography')
    if draw_coastline :
        ax.contour(x,y,c,
                   levels=[0.5],colors='black',linewidths=1.5,linestyles='dashed')

    # Longitude labels
    ax.set_xlim(imin-0.5,imax+0.5)
    ax.set_xticks(np.arange(imin,imax,line_skip))
    ax.set_xticks(np.arange(imin,imax),minor=True)
    ax.set_xticklabels(ax.get_xticks(), rotation = 45)

    # Latitude labels
    ax.set_ylim(jmin-0.5,jmax+0.5)
    ax.set_yticks(np.arange(jmin,jmax,line_skip))
    ax.set_yticks(np.arange(jmin,jmax),minor=True)

    for i in range(imin,imax,line_skip):
        ax.axvline(i+0.5,color='grey',linewidth=0.5)
    for j in range(jmin,jmax,line_skip):
        ax.axhline(j+0.5,color='grey',linewidth=0.5)

    cbar = fig.colorbar(pc,shrink=0.85)
    for j in range(jmin,jmax,label_skip):
        for i in range(imin,imax,label_skip):
                lab = f"({j},{i})"
                ax.text(i,j,lab,ha='center',va='center',fontsize=8)


"""
Do a flood fill of connected ocean regions taking into account
connectivity across cyclic boundaries including tripole
"""
def mask_flood(mask_in,seed):
    from skimage.morphology import flood
    
    # Create an extended grid with cyclic points included explicitly
    ny,nx = np.shape(mask_in.values)
    mask_ext = np.zeros((ny+1,nx+2))
    
    # Copy main field
    mask_ext[:-1,1:-1] = mask_in.values

    # E-W cyclic conditions except northern boundary
    mask_ext[:-1,0] = mask_in[:,-1].values
    mask_ext[:-1,-1] = mask_in[:,0].values

    # Northern tripole boundary
    for i in range(nx):
        mask_ext[-1,i+1] = mask_in[-1,nx-1-i]
    # These should always be zero since they are at the grid pole
    mask_ext[-1,0] = mask_ext[-1,-2]
    mask_ext[-1,nx+1] = mask_ext[-1,1]
    
    # Run the flood on the extended grid
    mask_ext = flood(mask_ext, seed, connectivity=1)
    
    # Subset to original mask size (drop cyclic points)
    da_flood=xr.full_like(mask_in,0)
    da_flood[:,:] = mask_ext[:-1,1:-1]
    
    da_flood = da_flood.assign_attrs({"fill_seed" : seed})
    
    return da_flood

def ice9(i, j, source, xcyclic=True, tripolar=True):
  """
  An iterative (stack based) implementation of "Ice 9".

  The flood fill starts at [j,i] and treats any positive value of "source" as
  passable. Zero and negative values block flooding.

  xcyclic = True allows cyclic behavior in the last index. (default)
  tripolar = True allows a fold across the top-most edge. (default)

  Returns an array of 0's and 1's.
  """
  wetMask = 0*source
  (nj,ni) = wetMask.shape
  stack = set()
  stack.add( (j,i) )
  while stack:
    (j,i) = stack.pop()
    if wetMask[j,i] or source[j,i] <= 0: continue
    wetMask[j,i] = 1
    if i>0: stack.add( (j,i-1) )
    elif xcyclic: stack.add( (j,ni-1) )
    if i<ni-1: stack.add( (j,i+1) )
    elif xcyclic: stack.add( (j,0) )
    if j>0: stack.add( (j-1,i) )
    if j<nj-1: stack.add( (j+1,i) )
    elif tripolar: stack.add( (j,ni-1-i) ) # Tri-polar fold
    return wetMask

def remove_isolated_points(mask_orig):

    # Create an array with corners as the third dimension
    mask_3D = mask_orig.expand_dims(dim={'neighbor':4}).copy() # SW
    mask_3D[1,:,:]=mask_orig.roll(lonh=-1).values              # SE
    mask_3D[2,:,:]=mask_orig.roll(lath=-1).values              # NW
    mask_3D[3,:,:]=mask_orig.roll(lonh=-1,lath=-1).values      # NE

    mask_nb = mask_orig[-1,::-1]
    mask_3D[2,-1,:] = mask_nb.values                           # Fix tripole boundary
    mask_3D[3,-1,:] = mask_nb.roll(lonh=-1).values
    
    # Find the minimum value of the corners
    mask_u = mask_3D.min(dim='neighbor').compute()
    
    # Now find the neighbors of the original points
    mask_3D[0,:,:]=mask_u.values                            # NE
    mask_3D[1,:,:]=mask_u.roll(lonh=1).values               # NW
    mask_3D[2,:,:]=mask_u.roll(lath=1).values               # SE
    mask_3D[3,:,:]=mask_u.roll(lonh=1,lath=1).values        # SW
    
    mask_3D[:,0,:] = 0.                                     # Fix first row to land
    
    # Find how many corners are true
    n = mask_3D.sum(dim='neighbor')
    
    # Return a mask for points that have at least one true neighbor    
    return mask_orig.where(n > 0,0)
 

def map_topo_latlon(df,lon_beg,lon_end,lat_beg,lat_end,zmax,zmin=0,
                    zvar='depth',mvar='mask'):

    # Find box using nominal lat/lon
    lonh = df.lonh
    lath = df.lath

    ibeg = np.abs(lonh - lon_beg).argmin().values
    jbeg = np.abs(lath - lat_beg).argmin().values
    iend = np.abs(lonh - lon_end).argmin().values
    jend = np.abs(lath - lat_end).argmin().values

    print('ibeg=',ibeg,' jbeg=',jbeg,
      'lat/lon beg=',df.geolat[jbeg,ibeg].values,df.geolon[jbeg,ibeg].values)
    print('iend=',iend,' jend=',jend,
      ' lat/lon end=',df.geolat[jend,iend].values,df.geolon[jend,iend].values)
    fig,ax=plt.subplots(figsize=(10,7.5),constrained_layout=True)
    plt.set_cmap('Blues')
    
    # extract subdomain
    lon2D = df.geolon.isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values
    lat2D = df.geolat.isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values
    dep2D = df[zvar].isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values
    wet2D = df[mvar].isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values

    elon2D = df.geolonb.shift(lonq=1).isel(lonq=slice(ibeg,iend),latq=slice(jbeg,jend)).values
    elat2D = df.geolatb.shift(latq=1).isel(lonq=slice(ibeg,iend),latq=slice(jbeg,jend)).values

    nlat=np.shape(lon2D)[0]
    nlon=np.shape(lon2D)[1]
#    print(nlon,nlat)

    pc=ax.pcolormesh(lon2D,lat2D,dep2D,vmin=zmin,vmax=zmax)

    for i in range(nlon) :
        ilab = '{0:d}'.format(i+ibeg)
        ax.text(lon2D[0,i],lat2D[0,i],ilab,color='red',ha='center',va='center')

    for j in range(nlat) :
        jlab = '{0:d}'.format(j+jbeg)
        ax.text(lon2D[j,0],lat2D[j,0],jlab,color='red',ha='center',va='center')
    

    for i in range(1,nlon) :
        for j in range(1,nlat) :
            dlab = '{0:.0f}'.format(dep2D[j,i])
            if ( wet2D[j,i] > 0)  :
                ax.text(lon2D[j,i],lat2D[j,i],dlab,ha='center',va='center')
            
    for i in range(nlon) :
        ax.plot(df.geolonb[jbeg:jend,i+ibeg],
                df.geolatb[jbeg:jend,i+ibeg],color='grey')
    for j in range(nlat) :
        ax.plot(df.geolonb[j+jbeg,ibeg:iend],
                df.geolatb[j+jbeg,ibeg:iend],color='grey')
        
    return ax

def inspect_topo(df,zvar,lon_beg,lon_end,lat_beg,lat_end,zmax,zmin=0,place=None):
    
    import numpy as np
    import xarray as xr
    import matplotlib.pyplot as plt
    from cartopy import crs as ccrs, feature as cfeature
    
    # Find box using nominal lat/lon
    lonh = df.lonh
    lath = df.lath

    ibeg,jbeg = mom6_point(df,lon_beg,lat_beg)
    iend,jend = mom6_point(df,lon_end,lat_end)
    
    # extract subdomain
    lon2D = df.geolon.isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values
    lat2D = df.geolat.isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values
    dep2D = df[zvar].isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values
    dep_mean2D = df['D_mean'].isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values
    dep_median2D = df['D_median'].isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values
    dep_min2D = df['D_min'].isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values
    dep_max2D = df['D_max'].isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values
    wet2D = df['mask'].isel(lonh=slice(ibeg,iend),lath=slice(jbeg,jend)).values

    elon2D = df.geolonb.shift(lonq=1).isel(lonq=slice(ibeg,iend),latq=slice(jbeg,jend)).values
    elat2D = df.geolatb.shift(latq=1).isel(lonq=slice(ibeg,iend),latq=slice(jbeg,jend)).values

    nlat=np.shape(lon2D)[0]
    nlon=np.shape(lon2D)[1]

    fig = plt.figure(figsize=(10,7.5),constrained_layout=True)
    ax = plt.axes(projection=ccrs.PlateCarree())
    

    plt.set_cmap('Blues')
    
#   plot the model topography as a color shaded map
    pc=ax.pcolormesh(lon2D,lat2D,dep2D,vmin=zmin,vmax=zmax,transform=ccrs.PlateCarree())

    # annotate with ij index labels
    for i in range(nlon) :
        ilab = '({0:d},{1:d})'.format(jbeg,i+ibeg)
        ax.text(lon2D[0,i],lat2D[0,i],ilab,
                color='red',ha='center',va='center',fontsize=6,
                transform=ccrs.PlateCarree())

    for j in range(1,nlat) :
        jlab = '({0:d},{1:d})'.format(j+jbeg,ibeg)
        ax.text(lon2D[j,0],lat2D[j,0],jlab,
                color='red',ha='center',va='center',fontsize=6,
                transform=ccrs.PlateCarree())
    
    # annotate with cell min/max and model depth
    for i in range(1,nlon) :
        for j in range(1,nlat) :
            
            if ( wet2D[j,i] > 0)  :
                dlab = '{0:.0f} \n $_{{{1:.0f}<{2:.0f}<{3:.0f}}}$'.format(dep2D[j,i],dep_min2D[j,i],dep_median2D[j,i],dep_max2D[j,i])
                ax.text(lon2D[j,i],lat2D[j,i],dlab,ha='center',va='center',fontsize=8,transform=ccrs.PlateCarree())                
            
    # annotate with grid lines
    for i in range(nlon) :
        ax.plot(df.geolonb[jbeg:jend,i+ibeg],
                df.geolatb[jbeg:jend,i+ibeg],color='blue',linewidth=0.5,
                transform=ccrs.PlateCarree())
    for j in range(nlat) :
        ax.plot(df.geolonb[j+jbeg,ibeg:iend],
                df.geolatb[j+jbeg,ibeg:iend],color='blue',linewidth=0.5,
                transform=ccrs.PlateCarree())
        
    if place != None :
        zlab = '{0:d}'.format(int(place['depth']))
        lat = place['lat']
        lon = place['lon']
        ax.scatter(lon,lat,marker='*',color='red')
        lon = lon + 0.01
        ax.text(lon,lat,zlab,color='red',ha='left',va='center',fontsize=8)
        
    ax.coastlines(color='grey')
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linestyle='dotted')
        
    return ax

def create_soc_topo_table():
    
    table = {
    'South Sandwich Trench' : {'lat' : -61.25, 'lon' : -23.75, 'depth' : 4600., 'width' : 70.},
    'Orkney Passage' :        {'lat' : -60.66, 'lon' : -40.75, 'depth' : 3200., 'width' : -1.},
    'Powell Basin' :          {'lat' : -60.53, 'lon' : -48.25, 'depth' : 2000., 'width' : -1.},
    'Shackleton F.Z.' :       {'lat' : -60.75, 'lon' : -56.75, 'depth' : 3500., 'width' : -1.},
    'Drake Passage' :         {'lat' : -57.50, 'lon' : -65.00, 'depth' : 3500., 'width' : -1.},
    'Georgia Passage' :       {'lat' : -56.00, 'lon' : -31.00, 'depth' : 3200., 'width' : -1.},
    'Shag Rocks Passage' :    {'lat' : -53.00, 'lon' : -49.00, 'depth' : 3200., 'width' : -1.},
    'Falkland Ridge Gap' :    {'lat' : -49.00, 'lon' : -35.00, 'depth' : 5100., 'width' : -1.},
    'Vema Sill' :             {'lat' : -28.66, 'lon' : -38.00, 'depth' : 4200., 'width' : -1.},
    'Hunter Channel' :        {'lat' : -35.00, 'lon' : -27.00, 'depth' : 4300., 'width' : -1.},
    'Meteor F.Z.' :           {'lat' : -35.50, 'lon' : -18.00, 'depth' : 3500., 'width' : -1.},
    'Walvis Passage' :        {'lat' : -36.00, 'lon' :  -7.00, 'depth' : 4200., 'width' : 50.},
    'Namib Col' :             {'lat' : -22.00, 'lon' :  -7.25, 'depth' : 3000., 'width' : 50.},
    'Cox F.Z.' :              {'lat' : -32.00, 'lon' : -12.15, 'depth' : 3600., 'width' : -1.},
    'Rio Grande F.Z.' :       {'lat' : -26.00, 'lon' : -14.75, 'depth' : 3900., 'width' : 13.},
    'Rio de Janeiro F.Z.' :   {'lat' : -22.50, 'lon' : -13.25, 'depth' : 3900., 'width' : 35.},
    'Bagration F.Z' :         {'lat' : -16.50, 'lon' : -13.50, 'depth' : 3800., 'width' : -1.},
    'Cardno F.Z.' :           {'lat' : -14.00, 'lon' : -13.50, 'depth' : 3500., 'width' : -1.},
    'Ascension F.Z.' :        {'lat' :  -8.00, 'lon' : -14.00, 'depth' : 3800., 'width' : 24.},
    'Guinea Rise' :           {'lat' :  -5.00, 'lon' :  -1.00, 'depth' : 4300., 'width' : -1.},
    'Chain F.Z.' :            {'lat' :  -1.00, 'lon' : -14.18, 'depth' : 4050., 'width' : 10.},
    'Romanche F.Z.' :         {'lat' :  -0.83, 'lon' : -13.75, 'depth' : 4350., 'width' : 10.},
    'St Pauls F.Z.' :         {'lat' :   0.25, 'lon' : -27.50, 'depth' : 3500., 'width' : -1.},
    'Four North F.Z.' :       {'lat' :   4.00, 'lon' : -27.50, 'depth' : 3500., 'width' : -1.},
    'Ceara Rise' :            {'lat' :   3.00, 'lon' : -39.42, 'depth' : 4300., 'width' : 300.},
    'Doldrums F.Z.' :         {'lat' :   8.23, 'lon' : -37.86, 'depth' : 4000., 'width' : -1.},
    'Vema F.Z.' :             {'lat' :  10.75, 'lon' : -40.92, 'depth' : 4650., 'width' : 20.},
    'Kane Gap' :              {'lat' :   9.66, 'lon' : -19.83, 'depth' : 4400., 'width' : -1.},
    '1520 F.Z.' :             {'lat' :  15.66, 'lon' : -48.83, 'depth' : 4400., 'width' : -1.},
    'Grenada Passage' :       {'lat' :  11.50, 'lon' : -62.00, 'depth' : 1000., 'width' : -1.},
    'St. Vincent Passage' :   {'lat' :  13.50, 'lon' : -61.00, 'depth' : 1000., 'width' : -1.},
    'St Lucia Passage' :      {'lat' :  14.25, 'lon' : -61.00, 'depth' : 1000., 'width' : -1.},
    'Dominica' :              {'lat' : -15.00, 'lon' : -61.25, 'depth' : 1085., 'width' : -1.},
    'Anegada Passage' :       {'lat' :  18.40, 'lon' : -64.20, 'depth' : 1800., 'width' : 15.},
    'Windward Passage' :      {'lat' :  19.00, 'lon' : -74.00, 'depth' : 1650., 'width' : 22.},
    'Yucatan Channel' :       {'lat' :  22.00, 'lon' : -86.0,  'depth' : 1650., 'width' : 22.},
    'Florida St.' :           {'lat' :  24.25, 'lon' : -80.50, 'depth' :  800., 'width' : 60.,},
    'Kane F.Z.' :             {'lat' :  23.67, 'lon' : -45.92, 'depth' : 3400., 'width' : -1.},
    'Atlantis F.Z.' :         {'lat' :  30.00, 'lon' : -41.92, 'depth' : 3200., 'width' : -1.},
    'Hayes F.Z.' :            {'lat' :  33.50, 'lon' : -38.50, 'depth' : 3300., 'width' : 14.},
    'Oceanographer F.Z.' :    {'lat' :  35.33, 'lon' : -36.50, 'depth' : 3100., 'width' : -1.},
    'Discovery Gap' :         {'lat' :  36.90, 'lon' : -16.66, 'depth' : 4675., 'width' : 50.},
    'St. of Gibralter' :      {'lat' :  35.92, 'lon' :  -5.75, 'depth' :  284., 'width' : 10.},
    'Gibbs F.Z.' :            {'lat' :  52.66, 'lon' : -35.50, 'depth' : 3675., 'width' : 10.},
    'Wyville-Thompson Ridge': {'lat' :  60.16, 'lon' :  -7.75, 'depth' :  600., 'width' : -1.},
    'Faroe Bank Channel' :    {'lat' :  61.50, 'lon' :  -8.50, 'depth' :  800., 'width' : 15.},
    'Iceland Faeroes Ridge' : {'lat' :  64.25, 'lon' : -12.00, 'depth' :  450., 'width' : -1.},
    'Denmark St.' :           {'lat' :  66.00, 'lon' : -28.00, 'depth' :  600., 'width' : -1.},
    'Tonga Trench' :          {'lat' : -10.00, 'lon' :-175.00, 'depth' : 5000., 'width' : -1.},
    'Samoa Passage' :         {'lat' :  -8.66, 'lon' :-169.00, 'depth' : 4800., 'width' : 200.},
    'Clipperton Passage' :    {'lat' :   0.00, 'lon' :-155.00, 'depth' : 4600., 'width' : 250.},
    'Vitiaz St.' :            {'lat' :  -5.80, 'lon' : 147.40, 'depth' : 1100., 'width' : 30.},
    'Mariana Ridge' :         {'lat' : -11.75, 'lon' :-139.00, 'depth' : 4450., 'width' : -1.},
    'Kyushu-Palau Ridge' :    {'lat' :  18.00, 'lon' : 135.00, 'depth' : 4600., 'width' : -1.},
    'Clarion Passage' :       {'lat' :  12.66, 'lon' :-165.66, 'depth' : 5200., 'width' : 40.},
    'Horizon Passage' :       {'lat' :  18.13, 'lon' :-169.18, 'depth' : 5000., 'width' : 10.},
    'Luzon St.' :             {'lat' :  20.00, 'lon' : 121.00, 'depth' : 2000., 'width' : -1.},
    'Wake Passage' :          {'lat' :  19.00, 'lon' : 169.00, 'depth' : 5250., 'width' : -1.},
    'Emporor Seamounts 1' :   {'lat' :  31.33, 'lon' : 174.00, 'depth' : 5200., 'width' : -1.},
    'Nintoku Passage' :       {'lat' :  39.00, 'lon' : 170.00, 'depth' : 4000., 'width' : -1.},
    'Shatsky Rise' :          {'lat' :  43.50, 'lon' : 169.00, 'depth' : 5100., 'width' : -1.},
    'Emporar Seamounts 2' :   {'lat' :  55.08, 'lon' : 164.58, 'depth' : 5370,  'width' : -1.},
    'Near St.' :              {'lat' :  53.00, 'lon' : 170.50, 'depth' : 2000., 'width' : -1.},
    'Amchitcka Pass' :        {'lat' :  50.50, 'lon' :-179.50, 'depth' : 1155., 'width' : -1.},
    'Amutka Pass' :           {'lat' :  51.00, 'lon' :-172.00, 'depth' :  430., 'width' : -1.},
    'Bering St.' :            {'lat' :  64.25, 'lon' :-171.60, 'depth' :   50., 'width' : 85.},
    'Chile Rise' :            {'lat' : -45.00, 'lon' : -80.00, 'depth' : 3600., 'width' :-1.},
    'Valdivia F.Z.' :         {'lat' : -41.00, 'lon' : -87.00, 'depth' : 3700., 'width' : -1.},
    'Nazca Ridge' :           {'lat' : -24.00, 'lon' : -84.00, 'depth' : 3900., 'width' : -1.},
    'Peru-Chile Trench' :     {'lat' : -18.00, 'lon' : -76.00, 'depth' : 4900., 'width' : -1.},
    'Gofar F.Z.' :            {'lat' :  -6.00, 'lon' :-105.00, 'depth' : 3500., 'width' : -1.},
    'Carnegie Ridge' :        {'lat' :  -1.00, 'lon' : -85.00, 'depth' : 2330., 'width' : 3.},
    'Equador Trench' :        {'lat' :  -0.35, 'lon' : -81.00, 'depth' : 2920., 'width' : 5.},
    'Cocos Ridge' :           {'lat' :   2.00, 'lon' : -90.00, 'depth':  2110., 'width' : -1.},
    'Sequieros F.Z.' :        {'lat' :   8.00, 'lon' :-102.00, 'depth' : 3200., 'width' : -1.},
    'Princess Elizabeth Trough':{'lat':-64.00, 'lon' :  83.00, 'depth' : 3500., 'width' : -1},
    'Austr-Antarc Discordance':{'lat': -48.00, 'lon' : 120.00, 'depth' : 4000., 'width' :-1.},
    'Crozet-Kerguelen Gap' :  {'lat' : -48.92, 'lon' :  57.83, 'depth' : 4400., 'width' :-1.},
    'Kerguelen-Amsterdam' :   {'lat' : -40.16, 'lon' :  77.42, 'depth' : 3200., 'width' : -1.},
    'Melville F.Z.' :         {'lat' : -37.00, 'lon' :  59.00, 'depth' : 4500., 'width' : -1.},
    'Atlantis II F.Z.' :      {'lat' : -31.75, 'lon' :  57.25, 'depth' : 5000., 'width' : -1.},
    'Mozambique Channel' :    {'lat' : -16.16, 'lon' :  41.25, 'depth' : 2500., 'width' : -1.},
    'Amirante Passage' :      {'lat' :  -9.50, 'lon' :  54.25, 'depth' : 4200., 'width' : -1.},
    'Ninetyeast Ridge 1' :    {'lat' : -10.00, 'lon' :  89.00, 'depth' : 3600., 'width' : -1.},
    'Ninetyeast Ridge 2' :    {'lat' :  -3.00, 'lon' :  89.20, 'depth' : 3500., 'width' : -1.},
    'Carlsberg Ridge' :       {'lat' : -10.50, 'lon' :  56.80, 'depth' : 4000., 'width' : -1.},
    'Bab El-Mandeb' :         {'lat' :  13.73, 'lon' :  42.50, 'depth' :  137., 'width' : 32.},
    'Lombok St.' :            {'lat' :  -8.90, 'lon' : 116.00, 'depth' :  350., 'width' : 22.},
    'Sumbas St' :             {'lat' :  -9.16, 'lon' : 120.16, 'depth' :  560., 'width' : -1.},
    'Sawu St.' :              {'lat' : -10.50, 'lon' : 121.10, 'depth' : 1000., 'width' : -1.},
    'Roti St;' :              {'lat' :  10.75, 'lon' : 122.50, 'depth' : 1000., 'width' : -1.},
    'Timor Sea' :             {'lat' :  -8.66, 'lon' : 130.40, 'depth' : 1400., 'width' : 80.},
    'Flores/Banda Sea' :      {'lat' :  -7.66, 'lon' : 122.58, 'depth' : 2500., 'width' : -1.},
    'Makassar St.' :          {'lat' :  -5.40, 'lon' : 118.00, 'depth' :  550., 'width' : -1.},
    'Lifamatola S.' :         {'lat' :  -1.83, 'lon' : 126.82, 'depth' : 2000., 'width' : -1.},
    'Halmahera basin' :       {'lat' :  -1.25, 'lon' : 129.00, 'depth' :  500., 'width' : -1.},
    'Batjan Basin' :          {'lat' :  -1.00, 'lon' : 126.75, 'depth' : 2300., 'width' : -1.},
    'Morotai Basin' :         {'lat' :  -2.83, 'lon' : 127.83, 'depth' : 2500., 'width' : -1.},
    'Celebes Sea' :           {'lat' :   5.00, 'lon' : 126.00, 'depth' : 1300., 'width' : -1.},
    'Jan Mayen F.Z.' :        {'lat' :  72.00, 'lon' :  -6.00, 'depth' : 1400., 'width' : -1.},
    'Mohns Ridge' :           {'lat' :  73.00, 'lon' :   7.00, 'depth' : 2500., 'width' : -1.},
    'Fram St.' :              {'lat' :  78.50, 'lon' :   5.00, 'depth' : 2700., 'width' : 450.},
    'Nansen Ridge' :          {'lat' :  86.00, 'lon' : 120.00, 'depth' : 1500., 'width' : -1.},
    'Alpha Mendeleyev Ridge' :{'lat' :  84.50, 'lon' :-152.00, 'depth' : 2600., 'width' : -1.}
    }

    return table
        
#SOC_Topo_table = create_soc_topo_table()

#print(SOC_Topo_table.keys())
#for n,key in enumerate(SOC_Topo_table.keys()) :
    print('%3d' % n + '%25s' % key + '%10.1f' % SOC_Topo_table[key]['depth'])
