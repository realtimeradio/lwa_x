import h5py
import sys
import numpy as np
import time

def get_bl_order(n_ants):
    """
    Return the order of baseline data output by a CASPER correlator
    X engine.
    
    Extracted from the corr package -- https://github.com/ska-sa/corr
    """
    order1, order2 = [], []
    for i in range(n_ants):
        for j in range(int(n_ants/2),-1,-1):
            k = (i-j) % n_ants
            if i >= k: order1.append((k, i))
            else: order2.append((i, k))
    order2 = [o for o in order2 if o not in order1]
    return tuple([o for o in order1 + order2])

def get_ant_names():
    """
    Generate a list of antenna names, where position
    in the list indicates numbering in the data files.
    """
    return ["foo"]*350

def create_header(h5):
    NANTS_DATA = 196
    NANTS = 350
    NTIMES = 10
    n_bls = (NANTS * (NANTS + 1)) / 2
    bls = np.array(get_bl_order(NANTS_DATA))

    header = h5.create_group("Header")
    header.create_dataset("Nants_data", dtype="<i8", data=NANTS_DATA)
    header.create_dataset("Nants_telescope", dtype="<i8", data=NANTS)
    header.create_dataset("Nbls",   dtype="<i8", data=n_bls)
    # Nblts needs populating by receiver
    #header.create_dataset("Nblts",  dtype="<i8", data=n_bls * NTIMES)
    header.create_dataset("Nfreqs", dtype="<i8", data=8192/4*3)
    header.create_dataset("Npols",  dtype="<i8", data=2)
    header.create_dataset("Nspws",  dtype="<i8", data=1)
    header.create_dataset("Ntimes", dtype="<i8", data=NTIMES)
    header.create_dataset("altitude",    dtype="<f8", data=0.0)
    header.create_dataset("ant_1_array", dtype="<i8", data=bls[:,0])
    header.create_dataset("ant_2_array", dtype="<i8", data=bls[:,1])
    header.create_dataset("antenna_diameters", dtype="<f8", data=14.0)
    header.create_dataset("antenna_names",     dtype="|S5", shape=(NANTS,), data=get_ant_names())
    header.create_dataset("antenna_numbers",   dtype="<i8", shape=(NANTS,), data=range(NANTS))
    header.create_dataset("antenna_positions", dtype="<f8", data=np.zeros([NANTS,3]))
    header.create_dataset("channel_width",     dtype="<f8", data=250e6 / 8192)
    header.create_dataset("freq_array",        dtype="<f8", shape=(1, 8192/4*3), data=np.linspace(0, 250e6, 8192/4*3))
    header.create_dataset("history",   data="%s: Template file created\n" % time.ctime())
    header.create_dataset("telescope", data="HERA")
    header.create_dataset("integration_time", dtype="<f8", data=10.0)
    header.create_dataset("latitude",    dtype="<f8", data=0.0)
    header.create_dataset("longitude",   dtype="<f8", data= 0.0)
    # lst_array needs populating by receiver. Should be center of integrations in radians
    #header.create_dataset("lst_array",   dtype="<f8", data=np.zeros(n_bls * NTIMES))
    header.create_dataset("object_name", data="zenith")
    header.create_dataset("phase_type",  data="drift")
    header.create_dataset("polarization_array", dtype="<i8", data=[-5, -6, -7, -8])
    header.create_dataset("spw_array",      dtype="<i8", data=0)
    header.create_dataset("telescope_name", data="HERA")
    # time_array needs populating by receiver (should be center of integrations in JD)
    #header.create_dataset("time_array", dtype="<f8", data=np.zeros(n_bls * NTIMES))
    # uvw_needs populating by receiver: uvw = xyz(ant2) - xyz(ant1). Units, metres.
    #header.create_dataset("uvw_array",  dtype="<f8", data=np.zeros([n_bls * NTIMES, 3]))
    header.create_dataset("vis_units",  data="uncalib")
    #header.create_dataset("zenith_dec", dtype="<f8", data=np.zeros(n_bls * NTIMES))
    #header.create_dataset("zenith_ra",  dtype="<f8", data=np.zeros(n_bls * NTIMES))
    # !Some! extra_keywords need to be computed for each file
    add_extra_keywords(header)


def add_extra_keywords(obj):
    extras = obj.create_group("extra_keywords")
    extras.create_dataset("cm_info", dtype="|S")
    extras.create_dataset("cmver", dtype="|S")
    extras.create_dataset("st_type", dtype="|S16")
    # The following keywords can only be set when the file is written
    extras.create_dataset("duration", dtype="<f8")
    extras.create_dataset("obs_id", dtype="<i8")
    extras.create_dataset("startt", dtype="<f8")
    extras.create_dataset("stopt",  dtype="<f8")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: make_hera_hdf5_template.py <filename>"
        exit()
    
    with h5py.File(sys.argv[1], 'w') as h5:
        create_header(h5)

    
