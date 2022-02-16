####
# title: ometiff2neuro.py
#
# language: python3
# date: 2022-02-16
# license: MIT
# author: viviana kwong, jason lu, dasun madhawa premathilaka, elmar bucher
#
# description:
#   script to render multi channel multi slice ometiff files into the neuroglancer software. 
#   + https://docs.openmicroscopy.org/ome-model/latest/
#   + https://github.com/google/neuroglancer
#########

# library
import argparse
import neuroglancer
import neuroglancer.cli
import tiffile


# functions
def ometiff2neuro(
        s_pathfile_tiff,
    ):
    '''
    docstring goes here.
    '''
    print(f'hello world: {s_pathfile_tiff}')

    # load tiff image as numpy array
    a_img =  tiffile.imread(s_pathfile_tiff)


# run the code
if __name__ == '__main__':
    o_parser = argparse.ArgumentParser(description='Script to render ome.tiff files into the neuroglancer software.')
    # request ometiff filename as command line argument
    o_parser.add_argument('ometiff', type=str, nargs=1, help='ome.tiff path/filename')
    # start neuroglancer
    neuroglancer.cli.add_server_arguments(o_parser)
    neuroglancer.cli.handle_server_arguments(o_parser.parse_args())
    viewer = neuroglancer.Viewer()
    with viewer.txn() as s:
        # render ometiff
        ometiff2neuro(o_parser.parse_args().ometiff[0])
    # print neuroglancer viewer url
    print(viewer)

