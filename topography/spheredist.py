def SphericalDistance(lonA, latA, lonB, latB):

'''
 Haversine formula (sin squared form)
   Calculate distance in meters between two points on idealized spherical earth
   given longitudes and lattudes in radians
'''

  #Angular separation 
  delta_lon = abs(lonA - lonB) % 2.0*np.pi  # map into 0 to 2pi
  delta_lat = latA - latB

  term1 = np.sin(0.5*delta_lat)**2
  term2 = cos(latA)*cos(latB)*sin(0.5*delta_lon)**2
  temp = term1 + term2
  gamma = 2.0*np.arctan2(np.sqrt(temp),np.sqrt(c1-temp))

  return gamma

