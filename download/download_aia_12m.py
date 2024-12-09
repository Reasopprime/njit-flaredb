''' a script to download the 12 minutes cadence SDO/AIA EUV/UV data corresponding to HMI data
'''

import os
import glob
import pandas as pd
import astropy.units as u
from sunpy.net import Fido
from sunpy.net import attrs as a
from datetime import datetime, timedelta

def download_aia_12m(hmi_dir, wave_length, jsoc_email, euv_cri = 12, uv_cri = 24):
    '''Return the save directory that stores AIA/EUV and UV data in 12m cadence corresponding to downloaded HMI data.
       euv_cri is the time window selected for AIA/EUV data, the defualt is 12 seconds.
       uv_cri is the time window selected for AIA/EUV data, the defualt is 24 seconds.
       The HMI sharp/sharp_cea data must be downloaded before applying this procedure.
    '''
    if os.path.exist(hmi_dir):
        hmi_files = glob.glob(os.path.join(hmi_dir,'*.magnetogram.fits'))
        hmi_files.sort()
        event_dir = os.path.abspath(os.path.join(hmi_dir, os.pardir))
        
        for aiawave in wave_length:

            if aiawave in [94, 131, 171, 193, 211, 304, 335]:
                aia_save_dir = os.path.join(event_dir, 'aia_euv', '12m', f'{str(aiawave)}')
                if os.path.exists(aia_save_dir):
                    print(f'AIA files exist, please check {aia_save_dir}')
                else:
                    os.makedirs(aia_save_dir, exist_ok=True)
                    for fn_hmi in hmi_files:
                        hmi_time = datetime.strptime(fn_hmi.split('TAI')[0][-16:-1], '%Y%m%d_%H%M%S')
                        query_sdo = []
                        count = 0
                        while len(query_sdo) == 0 and count <= int(euv_cri/12)*2-1:
                            start_time = hmi_time + (-1)**count*timedelta(seconds=int((count+1)/2)*12)
                            end_time = start_time
                            count += 1
                            query_sdo = Fido.search(a.Time(start_time, end_time),
                                                    a.Wavelength(int(aiawave)*u.angstrom),
                                                    a.jsoc.Series.aia_lev1_euv_12s,
                                                    a.jsoc.Notify(jsoc_email),
                                                    a.jsoc.Segment.image)
                        if len(query_sdo) == 1:
                            sdo_file = Fido.fetch(query_sdo, path = aia_save_dir)
                    print(f'{aia_save_dir} downloaded.')
            
            if aiawave in [1600, 1700]:
                aia_save_dir = os.path.join(event_dir, 'aia_uv', '12m', f'{str(aiawave)}')
                if os.path.exists(aia_save_dir):
                    print(f'AIA files exist, please check {aia_save_dir}')
                else:
                    os.makedirs(aia_save_dir, exist_ok=True)
                    for fn_hmi in hmi_files:
                        hmi_time = datetime.strptime(fn_hmi.split('TAI')[0][-16:-1], '%Y%m%d_%H%M%S')
                        query_sdo = []
                        count = 0
                        while len(query_sdo) == 0 and count <= int(uv_cri/24)*2-1:
                            start_time = hmi_time + (-1)**count*timedelta(seconds=int((count+1)/2)*24)
                            end_time = start_time
                            count += 1
                            query_sdo = Fido.search(a.Time(start_time, end_time),
                                                    a.Wavelength(int(aiawave)*u.angstrom),
                                                    a.jsoc.Series.aia_lev1_uv_24s,
                                                    a.jsoc.Notify(jsoc_email),
                                                    a.jsoc.Segment.image)
                        if len(query_sdo) == 1:
                            sdo_file = Fido.fetch(query_sdo, path = aia_save_dir)
                    print(f'{aia_save_dir} downloaded.')
