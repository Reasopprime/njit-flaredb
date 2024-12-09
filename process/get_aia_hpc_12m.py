'''a script to crop the full disk AIA image with same field of view of the corresponding HMI data in HPC coordinates
'''
import os
import glob
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy import ndimage
import sunpy.map
from copy import deepcopy
from reproject import reproject_exact
from astropy.wcs import WCS
import numpy as np
from sunpy.map.header_helper import make_heliographic_header
from sunpy.coordinates import RotatedSunFrame
from datetime import datetime, timedelta
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

def get_coord_hpc(file_dir):
    ''' get the coordinate of the center of the HPC HMI data
    '''
    sharp_map = sunpy.map.Map(file_dir)
    rota = np.deg2rad(sharp_map.meta['crota2'])
    xcen = sharp_map.meta['crval1'] + \
           sharp_map.meta['cdelt1']*np.cos(rota)*((sharp_map.meta['naxis1']+1)/2-sharp_map.meta['crpix1']) + \
           sharp_map.meta['cdelt2']*np.sin(rota)*((sharp_map.meta['naxis2']+1)/2-sharp_map.meta['crpix2'])
    ycen = sharp_map.meta['crval1'] + \
           sharp_map.meta['cdelt1']*np.sin(rota)*((sharp_map.meta['naxis1']+1)/2-sharp_map.meta['crpix1']) + \
           sharp_map.meta['cdelt2']*np.cos(rota)*((sharp_map.meta['naxis2']+1)/2-sharp_map.meta['crpix2'])
    hmit = datetime.strptime(file_dir.split('TAI')[0][-16:-1], '%Y%m%d_%H%M%S')
    del sharp_map
    return xcen*u.arcsec, ycen*u.arcsec, hmit

def get_aia_pixel_size(file_dir):
    ''' get the pixel size of the aia cropped image with the same field of view of HPC HMI data
    '''
    sharp_map = sunpy.map.Map(file_dir)
    hmi_res = 0.504040956
    aia_res = 0.600758016
    aia_width = np.round(sharp_map.meta['naxis1']*hmi_res/aia_res)
    aia_height = np.round(sharp_map.meta['naxis2']*hmi_res/aia_res)
    del sharp_map
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

def get_hpc_aia_12m(hmi_dir, aia_raw_dir, aia_save_dir, harpnumber):
    ''' hmi_dir: directory of the HMI HPC files 
    '''
    hmi_files = glob.glob(os.path.join(hmi_dir,'*.magnetogram.fits'))
    hmi_files.sort()

    if len(hmi_files)>0:
        aia_width_pixel, aia_height_pixel = get_aia_pixel_size(hmi_files[0])
      
        for fn_hmi in hmi_files:
            hmi_xcen, hmi_ycen, fn_hmi_time = get_coord_hpc(fn_hmi)
            hmi_map = sunpy.map.Map(fn_hmi)
            hmi_xycen_coord = SkyCoord(hmi_xcen, hmi_ycen, frame=hmi_map.coordinate_frame)
            fn_aia = find_closest_aia(fn_hmi_time, aia_raw_dir, t_cri = 60)
          
            if fn_aia:
              
                aia_map = sunpy.map.Map(fn_aia)
                aiat = datetime.strptime(fn_aia.split('Z')[0][-17:], '%Y-%m-%dT%H%M%S')
                durations = (aiat - fn_hmi_time).total_seconds()*u.second
                diffrot_point = SkyCoord(RotatedSunFrame(base=hmi_xycen_coord, duration=durations))
                aia_xycen = diffrot_point.transform_to(aia_map.coordinate_frame)
                aia_xycen_pixel = aia_map.wcs.world_to_pixel(aia_xycen)
                left_pixel = np.round(aia_xycen_pixel[0] - aia_width_pixel/2)
                bottom_pixel = np.round(aia_xycen_pixel[1] - aia_height_pixel/2)
                aia_submap = aia_map.submap([left_pixel, bottom_pixel]*u.pixel, width = aia_width_pixel*u.pixel, height = aia_height_pixel*u.pixel)
              
                fits_out_dir = os.path.join(aia_save_dir,'fits')
                os.makedirs(fits_out_dir, exist_ok=True)
                aia_sub_out_fits = os.path.join(fits_out_dir,os.path.basename(fn_aia).replace('.image_lev1.',f'.hpc_image.{harpnumber}.lev1.'))
                aia_submap.save(aia_sub_out_fits)
                del aia_map
                del aia_submap
        
