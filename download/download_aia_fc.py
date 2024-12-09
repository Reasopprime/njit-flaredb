''' a script to download the full cadence SDO/AIA EUV/UV data
'''

import os
import glob
import pandas as pd
import astropy.units as u
from sunpy.net import Fido
from sunpy.net import attrs as a
from datetime import datetime, timedelta

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

def download_aia_fc(flare_start_time, flare_end_time, jsoc_email, save_dir, wave_length):
    '''Return the save directory that stores AIA/EUV data in 12s cadence and UV data in 24s cadence
       Download time range start from 5 minutes before and 5 minutes after the flare start time
    '''
    if os.path.exists(save_dir):
        start_time = flare_start_time - timedelta(minutes=5)
        end_time = flare_end_time + timedelta(minutes=5)
        save_time = datetime.strftime(flare_start_time,'%Y-%m-%dT%H%M')

        for aiawave in wave_length:

            if aiawave in [94, 131, 171, 193, 211, 304, 335]:
                query_sdo = Fido.search(a.Time(start_time, end_time),
                                        a.Sample(12*u.second),
                                        a.Wavelength(int(aiawave)*u.angstrom),
                                        a.jsoc.Series.aia_lev1_euv_12s,
                                        a.jsoc.Notify(jsoc_email),
                                        a.jsoc.Segment.image)
                aia_save_dir = os.path.join(save_dir, f'{save_time}', 'aia_euv', '12s', f'{str(aiawave)}')
                os.makedirs(aia_save_dir, exist_ok=True)
                sdo_file = Fido.fetch(query_sdo, path = aia_save_dir)
                print(f'{sdo_file} downloaded.')

            if aiawave in [1600, 1700]:
                query_sdo = Fido.search(a.Time(start_time, end_time),
                                        a.Sample(24*u.second),
                                        a.Wavelength(int(aiawave)*u.angstrom),
                                        a.jsoc.Series.aia_lev1_uv_24s,
                                        a.jsoc.Notify(jsoc_email),
                                        a.jsoc.Segment.image)
                aia_save_dir = os.path.join(save_dir, f'{save_time}', 'aia_uv', '24s', f'{str(aiawave)}')
                os.makedirs(aia_save_dir, exist_ok=True)
                sdo_file = Fido.fetch(query_sdo, path = aia_save_dir)
                print(f'{sdo_file} downloaded.')
