####
# title: ometiff2neuro.py
#
# language: python3
# date: 2022-02-16
# license: MIT
# author: viviana kwong, jason lu, dasun madhawa premathilaka, elmar bucher
#
# installation:
#
# run:
#
#
# description:
#   script to render multi channel multi slice ometiff files into the neuroglancer software. 
#   + https://docs.openmicroscopy.org/ome-model/latest/
#   + https://github.com/google/neuroglancer
#########

# library
import argparse
from matplotlib import cm
import neuroglancer
import neuroglancer.cli
import numpy as np
from skimage import filters, io
import sys


# functions
def ometiff2neuro(
        o_state,
        s_pathfile_tiff,
        ls_channel_label = [
            'DNA1','PD1','TLR3','SOX10',
            'DNA2','CD163','CD3D','PDL1',
            'DNA3','CD4','ICOS','HLADPB1',
            'DNA4','CD8A','CD68','GZMB',
            'DNA5','CD40L','LAG3','HLAA',
            'DNA6','SQSTM','VIN','TIM3',
            'DNA7','LAMP1/CD107A','PDL1_2','PD1_2',
            'nuc_segement'
        ], # None gives hex label
        lli_channel_color = None,
        o_cm = cm.turbo,  # color map, used if no color list is None
        o_thresh = filters.threshold_li,  # None
        di_coor = {'c':1, 'z':0, 'y':2, 'x':3},
        di_nm = {'z':40, 'y':10, 'x':10},
        e_render = {0,'PD1','TLR3','CD3D','PDL1','CD8A','GZMB','VIN',28}   # layer integer or label, None renders all channels if the RAM can handle it.
    ):
    '''
    docstring goes here.
    '''
    print(f'\nprocessing: {s_pathfile_tiff}')

    # load tiff image as numpy array
    a_img =  io.imread(s_pathfile_tiff)
    print(f'image shape ({[m[0] for m in sorted(di_coor.items(), key=lambda n: n[1])]}): {a_img.shape}')

    # handle input
    i_channel = a_img.shape[di_coor['c']]  # channel count

    # handle e_render
    if e_render is None:
        e_render = set(range(i_channel))

    # handle channel labels
    if ls_channel_label is None:
        ls_channel_label = [hex(n) for n in range(i_channel)]
    elif len(ls_channel_label) != i_channel:
        sys.exit(f'Error @ ometiff2neuro : ls_channel_label shape (len(ls_channel_label)) does not match channel shape {i_channel}.')
    print(f'ls_channel_label: {ls_channel_label}')
        
    # handle channel color
    if lli_channel_color is None:
        lli_channel_color = []
        i_step = int(o_cm.N / i_channel)  # get step size to walk the color map
        for i_n in range(i_channel):
            lli_channel_color.append(o_cm(X=i_n*i_step, alpha=None, bytes=True))
    elif len(lli_channel_color) != i_channel:
        sys.exit(f'Error @ ometiff2neuro : lli_channel_color shape (len(ls_channel_label)) does not match channel shape {i_channel}.')
    print(f'lli_channel_color: {lli_channel_color}')

    # extract informiation and generate neuroglancer layer
    i_c = di_coor['c']
    for i_n in range(i_channel):
        s_label = ls_channel_label[i_n]
        print(f'check channel {i_n}/{i_channel}: {s_label}')
        if (i_n in e_render) or (s_label in e_render):
            print(f'rendering channel {i_n}/{i_channel}: {s_label}')

            # coordianet
            i_x = di_coor['x']
            i_y = di_coor['y']
            i_z = di_coor['z']
            if i_c == 0:
                a_channel = a_img[i_n,:,:,:]
            elif i_c == 1:
                a_channel = a_img[:,i_n,:,:]
                if di_coor['x'] > i_c: 
                    i_x -= 1
                if di_coor['y'] > i_c: 
                    i_y -= 1
                if di_coor['z'] > i_c: 
                    i_z -= 1
            elif i_c == 2:
                a_channel = a_img[:,:,i_n,:]
                if di_coor['x'] > i_c: 
                    i_x -= 1
                if di_coor['y'] > i_c: 
                    i_y -= 1
                if di_coor['z'] > i_c: 
                    i_z -= 1
            elif i_c == 3:
                a_channel = a_img[:,:,:,i_n]
            else:
                sys.exit(f'Error @ ometiff2neuro : the sorce code is broken. please fix the script.')

            # thresh
            if not (o_thresh is None):
                r_thresh = o_thresh(a_channel)
                a_channel[a_channel < r_thresh] = 0

            # segment
            ab_segment = a_channel > 0

            # color
            a_red = np.zeros(ab_segment.shape, dtype=np.uint8)
            a_red[ab_segment] = lli_channel_color[i_n][0]
            a_blue = np.zeros(ab_segment.shape, dtype=np.uint8)
            a_blue[ab_segment] = lli_channel_color[i_n][1]
            a_green = np.zeros(ab_segment.shape, dtype=np.uint8)
            a_green[ab_segment] = lli_channel_color[i_n][2]

            # result
            a_shape = np.array([a_red, a_blue, a_green])

            # generate neuroglancer object
            ls_name = ['c^', None, None, None]
            ls_name[i_x+1] = 'x'
            ls_name[i_y+1] = 'y'
            ls_name[i_z+1] = 'z'
            li_scale = [3, None, None, None]
            li_scale[i_x+1] = di_nm['x']
            li_scale[i_y+1] = di_nm['y']
            li_scale[i_z+1] = di_nm['z']
            state.layers.append(
                name = s_label,
                layer = neuroglancer.LocalVolume(
                    data = a_shape,
                    dimensions = neuroglancer.CoordinateSpace(
                        # rgb, x, y, z
                        names = ls_name,
                        scales = li_scale,
                        units = ['', 'nm', 'nm', 'nm'],
                    ),
                ),
                shader=
"""
void main() {
  emitRGB(
    vec3(toNormalized(getDataValue(0)),
    toNormalized(getDataValue(1)),
    toNormalized(getDataValue(2)))
  );
}
""",
            )


# run the code
if __name__ == '__main__':
    o_parser = argparse.ArgumentParser(description='Script to render ome.tiff files into the neuroglancer software.')
    # request path to ometiff and file name as command line argument
    o_parser.add_argument('ometiff', type=str, nargs=1, help='ome.tiff path/filename')
    # start neuroglancer
    neuroglancer.cli.add_server_arguments(o_parser)
    neuroglancer.cli.handle_server_arguments(o_parser.parse_args())
    viewer = neuroglancer.Viewer()
    with viewer.txn() as state:
        # render ometiff
        ometiff2neuro(
            o_state = state,
            s_pathfile_tiff = o_parser.parse_args().ometiff[0],
        )
    # print neuroglancer viewer url
    print(viewer)
