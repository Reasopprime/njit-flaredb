'''a script to make the quick look movie of AIA HPC/CEA data
'''
import os
import matplotlib.pyplot as plt
import sunpy.map
from sunpy.map.header_helper import make_heliographic_header
import warnings
from sunpy.map import MapSequence
import sunpy
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import matplotlib.colors as colors
import matplotlib.animation as animation
import sunpy.data.sample
import warnings
warnings.filterwarnings("ignore")


def MakeMovie_aia(data_dir, aia_wavelength, color_range, taitou, movie_name):
    '''taitou is the title of the output movie.
    '''
    sequence = sunpy.map.Map(os.path.join(data_dir, 'fits', '*.fits'), sequence=True)
    num_maps = len(sequence)
    aia_cmap = plt.get_cmap(f'sdoaia{aia_wavelength}')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection=sequence[0].wcs)  # WCSAxes for solar data
    im = sequence[0].plot(axes=ax, norm=colors.LogNorm(vmin=color_range[0], vmax=color_range[1]), cmap=aia_cmap)
    def update_frame(i):
        im.set_array(sequence[i].data)
        ax.set_title(f"{taitou}\n{sequence[i].date}")
    ani = animation.FuncAnimation(fig, update_frame, frames=num_maps, repeat=False)
    Writer = animation.writers['ffmpeg']   
    writer = Writer(fps=10, metadata=dict(artist='SunPy'), bitrate=1800)
    os.makedirs(os.path.join(data_dir, 'movie'), exist_ok=True)
    save_dir = os.path.join(data_dir, 'movie', f'{movie_name}')
    ani.save(save_dir, writer=writer)
