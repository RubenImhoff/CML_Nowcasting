# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:59:04 2019

@author: imhof002

Script to save the .dat-files, which contain the 15-min CML rainfall intensities
per grid cell location, as GeoTiff on the KNMI radar data grid.
"""

import h5py
import numpy as np
import pyproj
from osgeo import osr, ogr
from osgeo import gdal
import os
from scipy.spatial.kdtree import KDTree

############################
# Initial settings
############################

# Set the location of the h5-file
indir = 'f:\\CML\\NowcastingNL_ExtMask\\NowcastingNL\\2011\\RainMapsLinks15min'
#indir = 'm:\\My Documents\\PhD_internal\\CML_nowcasting\\subset_test'

# Set the output location + filename for the figure
outdir = 'f:\\CML\\NowcastingNL_ExtMask\\NowcastingNL\\RainMapsLinksGeoTiff15min\\2011'

# Set the filename of the grid with the interpolation grid 
coords_int_file = 'c:\\Documents\\PhD\\CML_nowcasting\\InterpolationGridextended.dat'

# Set the filename of the grid with the full radar (765 x 700) grid 
coords_2d_grid_file = 'c:\\Documents\\PhD\\CML_nowcasting\\radarcoordinaten_NL25_1km2_WGS84_full.dat'

# Filename h5 file --> We'll use one hdf5-file of the KNMI to convert the
# CML data into the same format
h5_filename = 'f:\\KNMI\\RAD_NL25_RAP_5min\\2012\\02\\RAD_NL25_RAP_5min_201202010045.h5'


#############################
# Set the function to plot the CML data on the 'normal' radar grid
#############################

def CML_to_radar_grid(filename, coords_int_file, coords_2d_grid_file, rows, cols):    
    # Make sure we don't get a scientific notation in the np arrays
    np.set_printoptions(suppress = True)
    
    # Open the 1D data on the interpolated grid (assumed that the location of each
    # cell corresponds to the cell in 'int_grid')
    f = np.genfromtxt(filename, delimiter = ' ', skip_header = 2, dtype = None)
    
    # Open the interpolation grid (first column is lon, second column is lat)
    int_grid = np.genfromtxt(coords_int_file, delimiter = ',', skip_header = 1, dtype = None)
    
    # Open the rectangular grid. This will be the final grid of the rainfall 
    # data for the nowcasting.
    radar_grid = np.genfromtxt(coords_2d_grid_file, delimiter = ',', skip_header = 1, dtype = None)
    
    # Make a tree of the radar dataset. We can query from this tree later on.
    tree = KDTree(list(zip(radar_grid[:,0], radar_grid[:,1])))
    
    # Find all locations in the radar_grid where int_grid has the same coordinates
    # as the radar grid.
    query_locs = tree.query(list(zip(int_grid[:,0], int_grid[:,1])))
    
    # Create the rectangular output grid on which we project the rainfall data.
    # We fill this with the KNMI nodata value (65535) for now.
    precip_intermediate = np.full(len(radar_grid), 65535)
    
    # Now, fill the precip_intermediate file with the values from f on the 
    # specified query locs.
    for i in range(0, len(int_grid)):
        precip_intermediate[query_locs[1][i]] = float(f[i])
        
    # Now reshape it to a 2D array
    precip = np.reshape(precip_intermediate, (rows,cols))

    return precip
    

############################
# Open the h5 file and get the necessary information
############################

# Open the file
ds_rainfall_file = h5py.File(h5_filename)
# Get the value out of the file
ds_rainfall_values = ds_rainfall_file['image1']['image_data'].value
# Set the value 65535 to NoData
ds_rainfall_values_v2 = np.where(ds_rainfall_values == 65535, np.NaN, ds_rainfall_values)
# Divide the values by 100, as they're saved as value * 100 [mm].
# Subsequently, multiply with 12 to get a rain rate in mm/h.
ds_rainfall_values_final = (ds_rainfall_values_v2 / 100.0) * 12.0

# Get the geographical information
geographic = ds_rainfall_file["geographic"]
# The projection
proj4str = geographic["map_projection"].attrs["projection_proj4_params"].decode()
pr = pyproj.Proj(proj4str)

# Get coordinates
latlon_corners = geographic.attrs["geo_product_corners"]
ll_lat = latlon_corners[1]
ll_lon = latlon_corners[0]
ur_lat = latlon_corners[5]
ur_lon = latlon_corners[4]
lr_lat = latlon_corners[7]
lr_lon = latlon_corners[6]
ul_lat = latlon_corners[3]
ul_lon = latlon_corners[2]

ll_x, ll_y = pr(ll_lon, ll_lat)
ur_x, ur_y = pr(ur_lon, ur_lat)
lr_x, lr_y = pr(lr_lon, lr_lat)
ul_x, ul_y = pr(ul_lon, ul_lat)
x1 = min(ll_x, ul_x)
y2 = min(ll_y, lr_y)
x2 = max(lr_x, ur_x)
y1 = max(ul_y, ur_y)


##########################
# Set the initial geodata (the data used for the transformation of the x- y-coordinates)
##########################
geodata = {}
geodata['x1'] = x1 
geodata['x2'] = x2 
geodata['y1'] = y1 
geodata['y2'] = y2 
geodata["xpixelsize"] = 1.0
geodata["ypixelsize"] = -1.0
geodata["projection"] = proj4str
geodata['lon1'] = min(ll_lon, ul_lon) 
geodata['lon2'] = max(lr_lon, ur_lon)
geodata['lat1'] = max(ul_lat, ur_lat)
geodata['lat2'] = min(ll_lat, lr_lat)
geodata["xpixelsize_wgs"] = (geodata['lon2'] - geodata['lon1'])/700
geodata["ypixelsize_wgs"] = (geodata['lat1'] - geodata['lat2'])/765
geodata["rows"] = geographic.attrs["geo_number_rows"][0]
geodata["cols"] = geographic.attrs["geo_number_columns"][0]

print("Number of rows is: ", geodata["rows"])
print("Number of cols is: ", geodata["cols"])

############################
# Create the GeoTIFF
############################

# Set a counter
counter = 0

# Start the loop
for infile in os.listdir(indir):
    print(infile)
    counter += 1
    print(float(counter)/float(len(os.listdir(indir))) * 100, ' % completed.')
    
    if infile.endswith('.dat'):
        # First get the rainfall data on a 2D grid with the function
        precip = CML_to_radar_grid(os.path.join(indir, infile), coords_int_file, coords_2d_grid_file, geodata['rows'], geodata['cols'])
    
        ##########################
        # Create the output raster
        ##########################
        
        # Initial settings
        dst_filename = os.path.join(outdir, (os.path.basename(infile)).split('_')[1].split('.')[0]+'.tif')
        x_pixels = int(geodata["cols"])
        y_pixels = int(geodata["rows"])
        pixel_size = geodata['xpixelsize']
        x_min = geodata['x1']
        y_max = geodata['y1']
        	
        driver = gdal.GetDriverByName('GTiff')
        
        # Put this in the gdal geoTIFF format
        dataset = driver.Create(
                dst_filename,
                x_pixels,
                y_pixels,
                1,
                gdal.GDT_Float32, )
        		
        dataset.SetGeoTransform((
                x_min,
                pixel_size,
                0,
                y_max,
                0,
                -pixel_size))
        
        # Set the projection
        projection = osr.SpatialReference()
        projection.ImportFromProj4(proj4str)
        
        # Save it
        dataset.SetProjection(projection.ExportToWkt())
        dataset.GetRasterBand(1).WriteArray(precip)
        dataset.FlushCache() # This is the way to write it to the disk
        
            # Open the file
        precip = None
        dst_filename = None

