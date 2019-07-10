import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse

parser = argparse.ArgumentParser(description='Read uvh5 files written out by the BDA catcher hashpipe pipeline and plot data',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('fname', type=str, help='uvh5 file to read')
args = parser.parse_args()

fp = h5py.File(args.fname,'r')
data = fp['Data']['visdata'].value


