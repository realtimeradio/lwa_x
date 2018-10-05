import h5py
import sys
import numpy as np
import time
import pickle

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
    return ["foo"]*352

def get_cm_info():
    from hera_mc import sys_handling
    h = sys_handling.Handling()
    return h.get_cminfo_correlator()

def get_antpos_enu(antpos, lat, lon, alt):
    import pyuvdata
    ecef = pyuvdata.ECEF_from_rotECEF(antpos, lon)
    enu  = pyuvdata.ENU_from_ECEF(ecef, lat, lon, alt)
    return enu

def create_header(h5, use_cm=False):
    """
    Create an HDF5 file with appropriate datasets in a "Header"
    data group.

    inputs: h5 -- an h5py File object
            use_cm -- boolean. If True, get valid data from the hera_cm
                      system. If False, just stuff the header with fake
                      data.
    """
    if use_cm:
        cminfo = get_cm_info()
        # add the enu co-ords
        cminfo["antenna_positions_enu"] = get_antpos_enu(cminfo["antenna_positions"], cminfo["cofa_lat"],
                                                             cminfo["cofa_lon"], cminfo["cofa_alt"])

    INSTRUMENT = "HERA"
    NANTS_DATA = 192
    NANTS = 352
    NCHANS = 2048 / 4. * 3
    ANT_DIAMETER = 14.0
    INT_TIME = 10.0
    n_bls = (NANTS * (NANTS + 1)) / 2
    bls = np.array(get_bl_order(NANTS_DATA))

    header = h5.create_group("Header")
    header.create_dataset("Nants_data", dtype="<i8", data=NANTS_DATA)
    header.create_dataset("Nants_telescope", dtype="<i8", data=NANTS)
    header.create_dataset("Nbls",   dtype="<i8", data=n_bls)
    header.create_dataset("Nblts",  dtype="<i8", data=0) # updated by receiver when file is closed
    header.create_dataset("Nfreqs", dtype="<i8", data=NCHANS)
    header.create_dataset("Npols",  dtype="<i8", data=2)
    header.create_dataset("Nspws",  dtype="<i8", data=1)
    header.create_dataset("Ntimes", dtype="<i8", data=0) # updated by receiver when file is closed
    header.create_dataset("ant_1_array", dtype="<i8", data=bls[:,0])
    header.create_dataset("ant_2_array", dtype="<i8", data=bls[:,1])
    header.create_dataset("antenna_diameters", dtype="<f8", data=ANT_DIAMETER)
    header.create_dataset("channel_width",     dtype="<f8", data=250e6 / (NCHANS / 3 * 4))
    header.create_dataset("freq_array",        dtype="<f8", shape=(1, NCHANS), data=np.linspace(0, 250e6, NCHANS))
    header.create_dataset("history",   data="%s: Template file created\n" % time.ctime())
    header.create_dataset("instrument", data=INSTRUMENT)
    header.create_dataset("integration_time", dtype="<f8", data=INT_TIME)
    header.create_dataset("object_name", data="zenith")
    header.create_dataset("phase_type",  data="drift")
    header.create_dataset("polarization_array", dtype="<i8", data=[-5, -6, -7, -8])
    header.create_dataset("spw_array",      dtype="<i8", data=0)
    header.create_dataset("telescope_name", data="HERA")
    header.create_dataset("vis_units",  data="uncalib")
    if use_cm:
        header.create_dataset("altitude",    dtype="<f8", data=cminfo['cofa_alt'])
        ant_pos = np.zeros([NANTS,3], dtype=np.int64)
        ant_pos_enu = np.zeros([NANTS,3], dtype=np.int64)
        ant_names = ["NONE"]*NANTS
        ant_nums = [-1]*NANTS
        for n, i in enumerate(cminfo["antenna_numbers"]):
            ant_pos[i]     = cminfo["antenna_positions"][n]
            ant_names[i]   = cminfo["antenna_names"][n].encode()
            ant_nums[i]    = cminfo["antenna_numbers"][n]
            ant_pos_enu[i] = cminfo["antenna_positions_enu"][n]
        header.create_dataset("antenna_names",     dtype="|S5", shape=(NANTS,), data=ant_names)
        header.create_dataset("antenna_numbers",   dtype="<i8", shape=(NANTS,), data=ant_nums)
        header.create_dataset("antenna_positions",   dtype="<i8", shape=(NANTS,3), data=ant_pos)
        header.create_dataset("antenna_positions_enu",   dtype="<i8", shape=(NANTS,3), data=ant_pos_enu)
        header.create_dataset("latitude",    dtype="<f8", data=cminfo["cofa_lat"])
        header.create_dataset("longitude",   dtype="<f8", data=cminfo["cofa_lon"])
    else:
        header.create_dataset("altitude",    dtype="<f8", data=0.0)
        header.create_dataset("antenna_names",     dtype="|S5", shape=(NANTS,), data=["NONE"]*NANTS)
        header.create_dataset("antenna_numbers",   dtype="<i8", shape=(NANTS,), data=range(NANTS))
        header.create_dataset("antenna_positions",   dtype="<i8", shape=(NANTS,3), data=np.zeros([NANTS,3]))
        header.create_dataset("antenna_positions_enu",   dtype="<i8", shape=(NANTS,3), data=np.zeros([NANTS,3]))
        header.create_dataset("latitude",    dtype="<f8", data=0.0)
        header.create_dataset("longitude",   dtype="<f8", data=0.0)

    # lst_array needs populating by receiver. Should be center of integrations in radians
    #header.create_dataset("lst_array",   dtype="<f8", data=np.zeros(n_bls * NTIMES))
    # time_array needs populating by receiver (should be center of integrations in JD)
    #header.create_dataset("time_array", dtype="<f8", data=np.zeros(n_bls * NTIMES))
    # uvw_needs populating by receiver: uvw = xyz(ant2) - xyz(ant1). Units, metres.
    #header.create_dataset("uvw_array",  dtype="<f8", data=np.zeros([n_bls * NTIMES, 3]))
    # !Some! extra_keywords need to be computed for each file
    if use_cm:
        add_extra_keywords(header, cminfo)
    else:
        add_extra_keywords(header)


def add_extra_keywords(obj, cminfo=None):
    extras = obj.create_group("extra_keywords")
    if cminfo is not None:
        extras.create_dataset("cmver", data=cminfo["cm_version"])
        extras.create_dataset("cminfo", data=pickle.dumps(cminfo))
    else:
        extras.create_dataset("cmver", data="generated-without-cminfo")
        extras.create_dataset("cminfo", data="generated-without-cminfo")
        
    extras.create_dataset("st_type", data="???")
    extras.create_dataset("duration", dtype="<f8", data=0.0) # filled in by receiver
    extras.create_dataset("obs_id", dtype="<i8", data=0)     # filled in by receiver
    extras.create_dataset("startt", dtype="<f8", data=0.0)   # filled in by receiver
    extras.create_dataset("stopt",  dtype="<f8", data=0.0)   # filled in by receiver

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Create a template HDF5 header file, optionally '\
                                     'using the correlator C+M system to get current meta-data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('output',type=str, help = 'Path to which the template header file should be output')
    parser.add_argument('-c', dest='use_cminfo', action='store_true', default=False,
                        help ='Use this flag to get up-to-date (hopefully) array meta-data from the C+M system')
    args = parser.parse_args()

    with h5py.File(args.output, "w") as h5:
        create_header(h5, args.use_cminfo)
