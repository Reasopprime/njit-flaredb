'''a script to project the full disk AIA image into cylindrical equal area (CEA) coordinates
   then crop into the same field of view of the corresponding HMI CEA data 
'''
import os
import psutil
import matplotlib.pyplot as plt
import glob
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy import ndimage
import sunpy.map
from copy import deepcopy
from reproject import reproject_exact
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import numpy as np
from sunpy.map.header_helper import make_heliographic_header
from sunpy.coordinates import RotatedSunFrame
from datetime import datetime, timedelta
import pandas as pd
import warnings
import gc
warnings.filterwarnings("ignore")

def get_coord_cea(file_dir):
    ''' get the coordinate of the center of the CEA HMI data
    '''
    cea_map = sunpy.map.Map(file_dir)
    cea_xcen_deg = cea_map.meta['crval1']*u.degree
    cea_ycen_deg = cea_map.meta['crval2']*u.degree
    hmit = datetime.strptime(file_dir.split('TAI')[0][-16:-1], '%Y%m%d_%H%M%S')
    del cea_map
    gc.collect()
    return cea_xcen_deg, cea_ycen_deg, hmit

def get_cea_pixel_size(file_dir):
    ''' get the pixel size of the aia cropped image with the same field of view of CEA HMI data
    '''
    cea_map = sunpy.map.Map(file_dir)
    aia_width = cea_map.meta['naxis1']
    aia_height = cea_map.meta['naxis2']
    del cea_map
    gc.collect()
    return aia_width, aia_height

def find_closest_aia(target_time, aia_dir, t_cri):
    ''' find the corresponding AIA 12m file to a target HPC HMI file
        t_cri is the searching range before and after the target HMI time, 
        usually t_cri is depending on the euv_cri and uv_cri in download_aia_12m.py
    '''
    aia_files = glob.glob(os.path.join(aia_dir, '*.image_lev1.fits'))
    aia_files.sort()
    aia_times = []
    for fn in aia_files:
        extract_time = datetime.strptime(fn.split('Z')[0][-17:], '%Y-%m-%dT%H%M%S')
        aia_times.append(extract_time)
    closest_index, closest_time = min(enumerate(aia_times), key=lambda x: abs(x[1] - target_time))
    time_diff = abs(closest_time - target_time)
    if time_diff.total_seconds() <= t_cri:
        return aia_files[closest_index]
    else:
        return False

def get_aia_cea_12m(hmi_dir, aia_raw_dir, aia_save_dir, harpnumber):
    ''' hmi_dir: directory of the HMI CEA files 
    '''
    hmi_files = glob.glob(os.path.join(hmi_dir,'*.magnetogram.fits'))
    hmi_files.sort()

    if len(hmi_files)>0:
        cea_width_pixel, cea_height_pixel = get_cea_pixel_size(hmi_files[0])
      
        for fn_hmi in hmi_files:
            cea_xcen, cea_ycen, fn_hmi_time = get_coord_cea(fn_hmi)
            #hmi_map = sunpy.map.Map(fn_hmi)
            #hmi_xycen_coord = SkyCoord(hmi_xcen, hmi_ycen, frame=hmi_map.coordinate_frame)
            fn_aia = find_closest_aia(fn_hmi_time, aia_raw_dir, t_cri = 60)
          
            if fn_aia:
              
                aia_map = sunpy.map.Map(fn_aia)
                shape = (6000,12000)
                carr_header = make_heliographic_header(aia_map.date,aia_map.observer_coordinate, shape, frame='carrington')
                aia_cea = aia_map.reproject_to(carr_header)
                del aia_map
                gc.collect

                xycen = SkyCoord(cea_xcen, cea_ycen, frame=aia_cea.coordinate_frame)
                xycen_pixel = aia_cea.wcs.world_to_pixel(xycen)
                left_pixel = np.round(xycen_pixel[0] - cea_width_pixel/2)
                bottom_pixel = np.round(xycen_pixel[1] - cea_height_pixel/2)
                aia_cea_sub = aia_cea.submap([left_pixel, bottom_pixel]*u.pixel, width = (cea_width_pixel-1)*u.pixel, height = (cea_height_pixel-1)*u.pixel)
                del aia_cea
                gc.collect

                aia_cea_sub_data = aia_cea_sub.data.astype(np.float32)
                aia_save = sunpy.map.Map(aia_cea_sub_data, aia_cea_sub.meta)
                del aia_cea_sub
                gc.collect

                fits_out_dir = os.path.join(aia_save_dir, 'fits')
                os.makedirs(fits_out_dir, exist_ok=True)
                aia_sub_out_fits = os.path.join(fits_out_dir,os.path.basename(fn_aia).replace('.image_lev1.',f'.cea_image.{harpnumber}.lev1.'))
                aia_save.save(aia_sub_out_fits)
                del aia_save
                gc.collect()
