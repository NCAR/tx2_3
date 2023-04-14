import numpy as np
import xarray as xr

import matplotlib.pyplot as plt

import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point

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
    
def map_mask(df,vmask,lonbeg,lonend,latbeg,latend,skip=10,
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
    if ( latend > 65. ) :
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
    
    for i in range(imin,imax):
        ax.plot(lone[:,i]+lo,late[:,i],
                color='grey',linewidth=0.5)
    for j in range(jmin,jmax):
        ax.plot(lone[j,:]+lo,late[j,:],
                color='grey',linewidth=0.5)

    cbar = fig.colorbar(pc,shrink=0.85)
    for j in range(jmin,jmax,skip):
        for i in range(imin,imax,skip):
            if ( latbeg < latc[j,i] and
                 latend > latc[j,i] and
                 lonbeg+lo < lonc[j,i]+lo and
                 lonend+lo > lonc[j,i]+lo ) :
                       
                lab = f"({j},{i})"
                ax.text(lonc[j,i]+lo,latc[j,i],
                        lab,ha='center',va='center',fontsize=8)


