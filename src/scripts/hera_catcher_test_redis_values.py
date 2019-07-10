import numpy as np
import redis
import time

r = redis.Redis("redishost")

hdf5template = '/home/deepthi/paper_gpu/src/hdf5header_template.txt'
configfile = '/home/deepthi/paper_gpu/src/scripts/bda_config_16ants_nobda.txt'
synctime = time.time()
inttime  = 64*2048
nfiles = 2
trigger = 1

pubchan = 'hashpipe://px1/0/set'

r.publish(pubchan, 'HDF5TPLT=%s' % hdf5template)
r.publish(pubchan, 'SYNCTIME=%d' % int(synctime))
r.publish(pubchan, 'INTTIME=%d' % inttime)
r.publish(pubchan, 'NFILES=%d' % nfiles)
r.publish(pubchan, 'DISKMING=9999');
r.publish(pubchan, 'BDACONFIG=%s' % configfile.split('/')[-1])

r.hset('hashpipe_bda','config', configfile)

baselines = {}
for n in range(4):
    baselines[n] = []

bdaconfig = np.loadtxt(configfile, dtype=np.int)
for i,t in enumerate(bdaconfig[:,2]):
    if (t==0): continue
    n = int(np.log(t)/np.log(2))
    if (n==4): n = 3
    baselines[n].append((bdaconfig[i,0], bdaconfig[i,1]))

for i in range(4):
    r.publish(pubchan, 'NBL%dSEC=%d'  % (2**(i+1), len(baselines[i])))

time.sleep(0.1)
r.publish(pubchan, 'TRIGGER=%d' % trigger)
