#########
# title: ometiff2neuro.py
#
# language: python3
# date: 2022-02-16
# license: MIT
# author: elmar bucherl, jason lu, viviana kwong
#
# installation:
#     conda create -n neuro python=3
#     conda activate neuro
#     pip install neuroglancer
#     pip install ipython  # for coding
#     pip install matplotlib  # for color maps
#     pip install scikit-image  # for loading images as numpy array and signal thresh
#     #pip install aicsimageio  # to extract ometiff metadata
#
# test dataset:
#     conda activate neuro
#     pip install synapseclient  # to download the data from synapse
#     synapse get -r syn26848775
#
# run:
#     python3 -i ometiff2neuro.py <path/filename.tiff>
#
# description:
#   script to render multi channel multi slice (ome)tiff files into the neuroglancer software.
#   + https://docs.openmicroscopy.org/ome-model/latest/
#   + https://github.com/google/neuroglancer
#########

# library
import argparse
from matplotlib import cm
import neuroglancer
import neuroglancer.cli
import numpy as np
from skimage import exposure, filters, io, util
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
        ], # None gives hex label.
        o_intensity_cm = cm.gray,  # color map, used for intensity. if None, no intensity layer.
        b_intensity_norm = True,
        o_thresh = filters.threshold_li,  # set None if segmetation mask data or already threshed values but only then!
        di_coor = {'c':1, 'z':0, 'y':2, 'x':3}, # dataset dependent. information would be in ometiff metadata.
        di_nm = {'z':200, 'y':108, 'x':108},  # microscope dependent. information would be in ometiff metadata.
        e_render = {0, 'CD4', 'CD8A', 'GZMB', 28},  # layer integer or label, None renders all channels if the RAM can handle it.
    ):
    '''
    input:
        o_state:
        s_pathfile_tiff:
        ls_channel_label:
        o_intensity_cm:
        b_intensity_norm:
        o_thresh:
        di_coor:
        di_nm:
        e_render:

    output:
        url

    description:
        docstring goes here.
    '''
    print(f'\nprocessing: {s_pathfile_tiff}')

    # load tiff image as numpy array
    a_img =  io.imread(s_pathfile_tiff)
    print(f'image shape ({[m[0] for m in sorted(di_coor.items(), key=lambda n: n[1])]}): {a_img.shape}')

    # channel count
    i_channel = a_img.shape[di_coor['c']]

    # handle channel labels
    if ls_channel_label is None:
        ls_channel_label = [hex(n) for n in range(i_channel)]
    elif len(ls_channel_label) != i_channel:
        sys.exit(f'Error @ ometiff2neuro : ls_channel_label shape (len(ls_channel_label)) does not match channel shape {i_channel}.')
    print(f'ls_channel_label: {ls_channel_label}')

    # handle render set
    if e_render is None:
        e_render = set(range(i_channel))

    # generate neuroglancer layer #
    i_c = di_coor['c']
    for i_n in range(i_channel):
        s_label = ls_channel_label[i_n]
        print(f'check channel {i_n}/{i_channel}: {s_label}')
        if (i_n in e_render) or (s_label in e_render):
            print(f'rendering channel {i_n}/{i_channel}: {s_label}')

            # extract data and coordinate columns
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
                sys.exit(f'Error @ ometiff2neuro : the source code is broken. please, fix the script.')


            # 3D rendering #
            # thresh data
            a_thresh = a_channel.copy()
            if not (o_thresh is None):
                r_thresh = o_thresh(a_thresh)
                a_thresh[a_thresh < r_thresh] = 0
            ab_thresh = a_thresh > 0

            # shape
            a_shape = np.zeros(ab_thresh.shape, dtype=np.uint32)
            a_shape[ab_thresh] = i_n + 1

            # generate neuroglancer object
            ls_name = [None, None, None]
            ls_name[i_x] = 'x'
            ls_name[i_y] = 'y'
            ls_name[i_z] = 'z'
            li_scale = [None, None, None]
            li_scale[i_x] = di_nm['x']
            li_scale[i_y] = di_nm['y']
            li_scale[i_z] = di_nm['z']
            state.layers.append(
                name = s_label,
                layer = neuroglancer.LocalVolume(
                    data = a_shape,
                    dimensions = neuroglancer.CoordinateSpace(
                        # rgb, x, y, z
                        names = ls_name,
                        scales = li_scale,
                        units = ['nm', 'nm', 'nm'],
                    ),
                ),
            )


            # expression intensity rendering #
            if not (o_intensity_cm is None):

                # normalize expression values by clip by two sigma and scale over the whole uint range
                if (b_intensity_norm):
                    i_min_clip = int(np.percentile(a_channel, 2.5))
                    i_max_clip = int(np.percentile(a_channel, 97.5))
                    a_clipped = np.clip(a_channel, a_min=i_min_clip, a_max=i_max_clip)
                    a_channel = exposure.rescale_intensity(a_clipped, in_range='image')  # 16 or 8[bit] normalized

                # translate intensity by color map
                a_8bit = util.img_as_ubyte(a_channel)
                a_intensity = o_intensity_cm(a_8bit, alpha=None, bytes=True)[:,:,:,0:3]

                # generate neuroglancer object
                ls_name = [None, None, None,'c^']
                ls_name[i_x] = 'x'
                ls_name[i_y] = 'y'
                ls_name[i_z] = 'z'
                li_scale = [None, None, None, 3]
                li_scale[i_x] = di_nm['x']
                li_scale[i_y] = di_nm['y']
                li_scale[i_z] = di_nm['z']
                state.layers.append(
                    name = s_label + '_intensity',
                    layer = neuroglancer.LocalVolume(
                        data = a_intensity,
                        dimensions = neuroglancer.CoordinateSpace(
                            # rgb, x, y, z
                            names = ls_name,
                            scales = li_scale,
                            units = ['nm', 'nm', 'nm', ''],
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
