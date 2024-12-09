'''a script to get the HMI data in HPC and CEA coordinates
   Crop the sharp data to a uniform size
   The segments of sharp and sharp_cea are change to Bx, By, Bz, continuum and magnetogram
   Bx and By are stored in the folder named as Bt (transverse field)
'''
import os
import glob
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import astropy.units as u
import shutil
import sunpy.map
from sunpy.map.header_helper import make_heliographic_header
from sunpy.coordinates import RotatedSunFrame
import warnings
warnings.filterwarnings("ignore")

def get_uni_size(sharp_dir):
    '''To get the uniform size of the sharp images with the center of the field of view remain fixed
    '''
    mag_files = glob.glob(os.path.join(sharp_dir,'*.magnetogram.fits'))
    mag_map = sunpy.map.Map(mag_files[0])
    x_min = len(mag_map.data[0])
    y_min = len(mag_map.data)
    del mag_map
    for i in range(1,len(mag_files)):
        mag_map = sunpy.map.Map(mag_files[i])
        x_size = len(mag_map.data[0])
        y_size = len(mag_map.data)
        x_min = x_min if x_min <= x_size else x_size
        y_min = y_min if y_min <= y_size else y_size
        del mag_map
    return [x_min, y_min]

def get_sharp_pro(sharp_raw_dir, sharp_pro_dir, uni_size):
    fn_mag = glob.glob(os.path.join(sharp_raw_dir, '*magnetogram.fits'))
    fn_mag.sort()
    os.makedirs(os.path.join(sharp_pro_dir,'sharp'), exist_ok=True)
    os.makedirs(os.path.join(sharp_pro_dir,'sharp','magnetogram','fits'), exist_ok=True)
    os.makedirs(os.path.join(sharp_pro_dir,'sharp','continuum','fits'), exist_ok=True)
    os.makedirs(os.path.join(sharp_pro_dir,'sharp','Bz','fits'), exist_ok=True)
    os.makedirs(os.path.join(sharp_pro_dir,'sharp','Bt','fits'), exist_ok=True)
    for fn_num in range(0,len(fn_mag)):
        mag_map = sunpy.map.Map(fn_mag[fn_num])
        con_map = sunpy.map.Map(fn_mag[fn_num].replace("magnetogram","continuum"))
        fld_map = sunpy.map.Map(fn_mag[fn_num].replace("magnetogram","field"))
        inc_map = sunpy.map.Map(fn_mag[fn_num].replace("magnetogram","inclination"))
        azi_map = sunpy.map.Map(fn_mag[fn_num].replace("magnetogram","azimuth"))
        
        mag_rot_map = sunpy.map.Map(np.rot90(mag_map.data,2), mag_map.meta)
        left_pixel = np.round(len(mag_map.data[0])/2 - uni_size[0]/2)
        bottom_pixel = np.round(len(mag_map.data)/2 - uni_size[1]/2)
        mag_sub = mag_rot_map.submap([left_pixel, bottom_pixel]*u.pixel, width = (uni_size[0]-1)*u.pixel, height = (uni_size[1]-1)*u.pixel)
        mag_fits_out = os.path.join(sharp_pro_dir,'sharp','magnetogram','fits',os.path.basename(fn_mag[fn_num]))
        mag_sub = sunpy.map.Map(mag_sub.data.astype('float32'), mag_sub.meta)
        mag_sub.save(mag_fits_out)
        del mag_rot_map, mag_sub

        con_rot_map = sunpy.map.Map(np.rot90(con_map.data,2), con_map.meta)
        con_sub = con_rot_map.submap([left_pixel, bottom_pixel]*u.pixel, width = (uni_size[0]-1)*u.pixel, height = (uni_size[1]-1)*u.pixel)
        con_fits_out = os.path.join(sharp_pro_dir,'sharp','continuum','fits',os.path.basename(fn_mag[fn_num].replace("magnetogram","continuum")))
        con_sub = sunpy.map.Map(con_sub.data.astype('float32'), con_sub.meta)
        con_sub.save(con_fits_out)
        del con_map, con_rot_map, con_sub
        
        bx = fld_map.data * np.sin(inc_map.data * u.degree) * np.sin(azi_map.data * u.degree)
        by = - fld_map.data * np.sin(inc_map.data * u.degree) * np.cos(azi_map.data * u.degree)
        bz = fld_map.data * np.cos(inc_map.data * u.degree)
        del fld_map, inc_map, azi_map

        bx_map = sunpy.map.Map(np.rot90(bx,2), mag_map.meta)
        bx_sub = bx_map.submap([left_pixel, bottom_pixel]*u.pixel, width = (uni_size[0]-1)*u.pixel, height = (uni_size[1]-1)*u.pixel)
        bx_fits_out = os.path.join(sharp_pro_dir,'sharp','Bt','fits',os.path.basename(fn_mag[fn_num].replace("magnetogram","Bx")))
        bx_sub = sunpy.map.Map(bx_sub.data.astype('float32'), bx_sub.meta)
        bx_sub.save(bx_fits_out)
        del bx_map, bx_sub

        by_map = sunpy.map.Map(np.rot90(by,2), mag_map.meta)
        by_sub = by_map.submap([left_pixel, bottom_pixel]*u.pixel, width = (uni_size[0]-1)*u.pixel, height = (uni_size[1]-1)*u.pixel)
        by_fits_out = os.path.join(sharp_pro_dir,'sharp','Bt','fits',os.path.basename(fn_mag[fn_num].replace("magnetogram","By")))
        by_sub = sunpy.map.Map(by_sub.data.astype('float32'), by_sub.meta)
        by_sub.save(by_fits_out)
        del by_map, by_sub

        bz_map = sunpy.map.Map(np.rot90(bz,2), mag_map.meta)
        bz_sub = bz_map.submap([left_pixel, bottom_pixel]*u.pixel, width = (uni_size[0]-1)*u.pixel, height = (uni_size[1]-1)*u.pixel)
        bz_fits_out = os.path.join(sharp_pro_dir,'sharp','Bz','fits',os.path.basename(fn_mag[fn_num].replace("magnetogram","Bz")))
        bz_sub = sunpy.map.Map(bz_sub.data.astype('float32'), bz_sub.meta)
        bz_sub.save(bz_fits_out)
        del bz_map, bz_sub, mag_map

def get_cea_pro(cea_raw_dir, cea_pro_dir):
    fn_mag = glob.glob(os.path.join(cea_raw_dir, '*magnetogram.fits'))
    fn_mag.sort()
    os.makedirs(os.path.join(cea_pro_dir,'cea'), exist_ok=True)
    os.makedirs(os.path.join(cea_pro_dir,'cea','magnetogram','fits'), exist_ok=True)
    os.makedirs(os.path.join(cea_pro_dir,'cea','continuum','fits'), exist_ok=True)
    os.makedirs(os.path.join(cea_pro_dir,'cea','Bz','fits'), exist_ok=True)
    os.makedirs(os.path.join(cea_pro_dir,'cea','Bt','fits'), exist_ok=True)
    for fn_num in range(0,len(fn_mag)):

        mag_fits_out = os.path.join(cea_pro_dir,'cea','magnetogram','fits',os.path.basename(fn_mag[fn_num]))
        shutil.copy(fn_mag[fn_num], mag_fits_out)

        con_fits_out = os.path.join(cea_pro_dir,'cea','continuum','fits',os.path.basename(fn_mag[fn_num].replace("magnetogram","continuum")))
        shutil.copy(fn_mag[fn_num].replace("magnetogram","continuum"), con_fits_out)

        bx_fits_out = os.path.join(cea_pro_dir,'cea','Bt','fits',os.path.basename(fn_mag[fn_num].replace("magnetogram","Bx")))
        bx_map = sunpy.map.Map(fn_mag[fn_num].replace("magnetogram","Bp"))
        bx_map = sunpy.map.Map(bx_map.data.astype('float32'), bx_map.meta)
        bx_map.save(bx_fits_out)

        bt_map = sunpy.map.Map(fn_mag[fn_num].replace("magnetogram","Bt"))
        by_data = - bt_map.data
        by_map = sunpy.map.Map(by_data.astype('float32'), bt_map.meta)
        by_fits_out = os.path.join(cea_pro_dir,'cea','Bt','fits',os.path.basename(fn_mag[fn_num].replace("magnetogram","By")))
        by_map.save(by_fits_out)

        bz_fits_out = os.path.join(cea_pro_dir,'cea','Bz','fits',os.path.basename(fn_mag[fn_num].replace("magnetogram","Bz")))
        bz_map = sunpy.map.Map(fn_mag[fn_num].replace("magnetogram","Br"))
        bz_map = sunpy.map.Map(bz_map.data.astype('float32'), bz_map.meta)
        bz_map.save(bz_fits_out)
        del bt_map, bx_map, by_map, bz_map
