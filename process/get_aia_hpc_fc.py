'''a script to crop the full disk AIA image in full cadence
   with same field of view of the corresponding HMI data in HPC coordinates
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

def get_flare_start_end_time(flare_list_file, flare_index):
    '''Return the flare start and end time of an event in the significant flare list
       Flare index start for 1 to 103
    '''
    flare_event_list = pd.read_csv(flare_list_file, header=None)
    start_time = flare_event_list.iloc[flare_index -1, 2][:16]
    flare_start_time = datetime.strptime(f'{start_time}', '%Y-%m-%d %H:%M')
    end_time = flare_event_list.iloc[flare_index -1, 4][:16]
    flare_end_time = datetime.strptime(f'{end_time}', '%Y-%m-%d %H:%M')
    return flare_start_time, flare_end_time

def get_hmi_target(flare_start_time, hmi_dir):
    '''find the most recent HMI file to the first AIA full cadence image
    '''
    hmi_start = flare_start_time - timedelta(minutes=flare_start_time.minute % 12)
    hmi_file = glob.glob(os.path.join(hmi_dir,f'*{hmi_start.strftime("%Y%m%d_%H%M%S")}_TAI.magnetogram.fits'))
    while not hmi_file:
        hmi_start = hmi_start - timedelta(minutes=12)
        hmi_file = glob.glob(os.path.join(hmi_dir,f'*{hmi_start.strftime("%Y%m%d_%H%M%S")}_TAI.magnetogram.fits'))
    return hmi_file, hmi_start

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
    del sharp_map
    return xcen*u.arcsec, ycen*u.arcsec

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

def get_hpc_aia_12s(hmi_dir, aia_raw_dir, aia_save_dir, flare_start_time, harpnumber):
  
    hmi_files = glob.glob(os.path.join(hmi_dir,'*.magnetogram.fits'))
    hmi_files.sort()

    if len(hmi_files)>0:
        hmi_target_file, hmi_target_time = get_hmi_target(flare_start_time, hmi_dir)
        aia_width_pixel, aia_height_pixel = get_aia_pixel_size(hmi_target_file)
        target_xcen, target_ycen = get_coord_hpc(hmi_target_file)
        hmi_target_map = sunpy.map.Map(hmi_target_file)
        point0 = SkyCoord(target_xcen, target_ycen, frame=hmi_target_map.coordinate_frame)
      
        aia_files = glob.glob(os.path.join(aia_raw_dir,'*.image_lev1.fits'))
        aia_files.sort()
        for fn_aia in aia_files:
            aia_map = sunpy.map.Map(fn_aia)
            aiat = datetime.strptime(fn_aia.split('Z')[0][-17:], '%Y-%m-%dT%H%M%S')
            durations = (aiat - hmi_target_time).total_seconds()*u.second
            diffrot_point = SkyCoord(RotatedSunFrame(base=point0, duration=durations))
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
