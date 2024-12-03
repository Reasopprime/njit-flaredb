''' a script to download the SDO/HMI sharp/sharp_cea LOS/vector magnetograms and continuum
'''

import os
import glob
import pandas as pd
import astropy.units as u
from sunpy.net import Fido
from sunpy.net import attrs as a
from datetime import datetime, timedelta

def get_flare_start_time_harpnum(flare_list_file, flare_index):
    '''Return the flare start and end time of an event in the significant flare list
       Flare index start for 1 to 103
    '''
    flare_event_list = pd.read_csv(flare_list_file, header=None)
    start_time = flare_event_list.iloc[flare_index -1, 2][:16]
    flare_start_time = datetime.strptime(f'{start_time}', '%Y-%m-%d %H:%M')
    harpnum = str(flare_event_list.iloc[flare_index -1, 5]).strip()
    return flare_start_time, harpnum

def download_sharp(flare_start_time, harp_number, jsoc_email, save_dir):
    '''Return the save directory that stores 32 hours of sharp data
       Format of the input flare_start_time: 'YYYY-MM-DD HH:MM'
       Download time range start from 24 hours before and 8 hours after the flare start time
    '''
    if os.path.exists(save_dir):
        start_time = flare_start_time - timedelta(minutes=flare_start_time.minute % 12) - timedelta(days=1)
        end_time = start_time + timedelta(hours=32)
        save_time = datetime.strftime(flare_start_time,'%Y-%m-%dT%H%M')

        result_sharp = Fido.search(a.Time(start_time, end_time),
                                   a.Sample(12*u.minute),
                                   a.jsoc.Series("hmi.sharp_720s"),
                                   a.jsoc.PrimeKey("HARPNUM", harp_number),
                                   a.jsoc.Notify(jsoc_email),
                                   a.jsoc.Segment("field") & 
                                   a.jsoc.Segment('inclination') &
                                   a.jsoc.Segment('azimuth') &
                                   #a.jsoc.Segment('disambig') &         # disambig solved since 15 Jan 2014
                                   a.jsoc.Segment('magnetogram') &
                                   a.jsoc.Segment('continuum'))
        sharp_save_dir = os.path.join(save_dir, 'sharp', f'{save_time}')
        os.makedirs(sharp_save_dir, exist_ok=True)
        file_sharp = Fido.fetch(result_sharp, path=sharp_save_dir, max_conn=10)
        
return file_sharp

def download_cea(flare_start_time, harp_number, jsoc_email, save_dir):
    '''Return the save directory of the file that stores 32 hours of sharp_cea data
       Format of the input flare_start_time: 'YYYY-MM-DD HH:MM'
       Download time range start from 24 hours before and 8 hours after the flare start time
    '''
    if os.path.exists(save_dir):
        start_time = flare_start_time - timedelta(minutes=flare_start_time.minute % 12) - timedelta(days=1)
        end_time = start_time + timedelta(hours=32)
        save_time = datetime.strftime(flare_start_time,'%Y-%m-%dT%H%M')

        result_cea = Fido.search(a.Time(start_time, end_time),
                                 a.Sample(12*u.minute),
                                 a.jsoc.Series("hmi.sharp_cea_720s"),
                                 a.jsoc.PrimeKey("HARPNUM", harp_number),
                                 a.jsoc.Notify(jsoc_email),
                                 a.jsoc.Segment("Br") & 
                                 a.jsoc.Segment('Bt') &
                                 a.jsoc.Segment('Bp') &
                                 a.jsoc.Segment('magnetogram') &
                                 a.jsoc.Segment('continuum'))
        cea_save_dir = os.path.join(save_dir, 'sharp_cea', f'{save_time}')
        os.makedirs(cea_save_dir, exist_ok=True)
        file_cea = Fido.fetch(result_cea, path=cea_save_dir, max_conn=10)

return file_cea
