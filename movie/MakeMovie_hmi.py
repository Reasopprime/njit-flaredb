'''a script to make the quick look movie of HMI HPC/CEA data
'''

import os
import matplotlib.pyplot as plt
from scipy import ndimage
import sunpy.map
import numpy as np
import warnings
from sunpy.map import MapSequence

import sunpy
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import matplotlib.colors as colors
import matplotlib.animation as animation
from astropy.time import Time
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.data.sample
import warnings
warnings.filterwarnings("ignore")

def MakeMovie_hmi(data_dir, if_Bt, if_continuum, color_range, taitou, movie_name):
    '''if_Bt: set this value as Bt to make transverse field movie using Bx and By
       if_continuum: set this value as continuum to make continuum movie
       taitou: title of the output movie
    '''
    if if_Bt == 'Bt':
        sequence_bx = sunpy.map.Map(os.path.join(data_dir, 'fits', '*Bx.fits'), sequence=True)
        sequence_by = sunpy.map.Map(os.path.join(data_dir, 'fits', '*By.fits'), sequence=True)
        bt_maps = []
        for bx_map, by_map in zip(sequence_bx, sequence_by):
            bt_data = np.sqrt(bx_map.data**2 + by_map.data**2)
            bt_map = sunpy.map.Map(bt_data, bx_map.meta)
            bt_maps.append(bt_map)
            sequence = MapSequence(bt_maps)
    else:
        sequence = sunpy.map.Map(os.path.join(data_dir, 'fits', '*.fits'), sequence=True)
    
    num_maps = len(sequence)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=sequence[0].wcs) 
    if if_continuum == 'continuum':
        color_range[0] = np.min(sequence[0].data)
        color_range[1] = np.max(sequence[0].data)
    im = sequence[0].plot(axes=ax, norm=ImageNormalize(vmin=color_range[0],vmax=color_range[1]))
    def update_frame(i):
        im.set_array(sequence[i].data)
        ax.set_title(f"{taitou}\n{sequence[i].date}")
    ani = animation.FuncAnimation(fig, update_frame, frames=num_maps, repeat=False)
    Writer = animation.writers['ffmpeg']   
    writer = Writer(fps=10, metadata=dict(artist='SunPy'), bitrate=1800)
    os.makedirs(os.path.join(data_dir, 'movie'), exist_ok=True)
    save_dir = os.path.join(data_dir, 'movie', f'{movie_name}')
    ani.save(save_dir, writer=writer)
