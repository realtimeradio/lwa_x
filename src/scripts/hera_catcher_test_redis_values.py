import numpy as np
import redis
import time

r = redis.Redis("redishost")

hdf5template = '/tmp/template.h5'
configfile = '/tmp/config.txt'
synctime = time.time()
inttime  = 64*2048
nfiles = 2
trigger = 1

pubchan = 'hashpipe://hera-sn1.corr.hera.pvt/0/set'

r.publish(pubchan, 'HDF5TPLT=%s' % hdf5template)
r.publish(pubchan, 'SYNCTIME=%d' % int(synctime))
r.publish(pubchan, 'INTTIME=%d' % inttime)
r.publish(pubchan, 'NFILES=%d' % nfiles)
r.publish(pubchan, 'DISKMING=9999');
r.publish(pubchan, 'BDACONF=%s' % configfile)

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

for v in ['NETWAT', 'NETREC', 'NETPRC']:
    r.publish(pubchan, '%sMN=99999' % (v))
    r.publish(pubchan, '%sMX=0' % (v))

r.publish(pubchan, 'NETDRPTL=0')

# Release nethread hold
r.publish(pubchan, 'CNETHOLD=0')
