#!/usr/bin/env python

import redis
import time
import argparse
import subprocess
import numpy as np

perf_tweaker = 'tweak-perf-sn.sh'
init = 'hera_catcher_init.sh'
python_source_cmd = ['source', '~/hera-venv/bin/activate']
template_cmd = ['hera_make_hdf5_template.py']
config = '/home/hera/dgorthi/paper_gpu/src/bdaconfig/test_bda_192ants_nobda.txt'

def run_on_hosts(hosts, cmd, user=None, wait=True):
    if isinstance(cmd, str):
        cmd = [cmd]
    p = []
    for host in hosts:
        if user is None:
            p += [subprocess.Popen(['ssh', '%s' % (host)] + cmd)]
        else:
            p += [subprocess.Popen(['ssh', '%s@%s' % (user, host)] + cmd)]
    if wait:
        for pn in p:
            pn.wait()

parser = argparse.ArgumentParser(description='Initialize the catcher machine and start hashpipe',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('host', type=str, help='Host to intialize')
parser.add_argument('-r', dest='redishost', type=str, default='redishost', help='Host serving redis database')
parser.add_argument('-t', dest='hdf5template', type=str, default='/tmp/template.h5', 
                    help='Place to put HDF5 header template file')
parser.add_argument('-c', dest='bdaconfig', type=str, default='/tmp/config.h5',
                    help='Place to put BDA config file')
parser.add_argument('--runtweak', dest='runtweak', action='store_true', default=False,
                    help='Run the tweaking script %s on X-hosts prior to starting the correlator' % perf_tweaker)

args = parser.parse_args()

r = redis.Redis(args.redishost)

## Run performance tweaking script
#if args.runtweak:
#    run_on_hosts([args.host], perf_tweaker, user='root', wait=True)
#
# Start Catcher
print ('Start catcher hashpipe script')
run_on_hosts([args.host], ['cd', '/data;', init, '0'], wait=False)

# Start hashpipe<->redis gateways
print('Establish redis gateway')
cpu_mask = '0x0004'
run_on_hosts([args.host], ['taskset', cpu_mask, 'hashpipe_redis_gateway.rb', '-g', args.host, '-i', '0'])

# Wait for the gateways to come up
time.sleep(15)

# Generate the meta-data template
print('Create template file')
run_on_hosts([args.host], python_source_cmd + [';'] + template_cmd + [args.hdf5template], wait=True)

# Copy config file to location request
print('Copy config file')
run_on_hosts([args.host], ['cp', config, args.bdaconfig], wait=True)

#Configure runtime parameters
catcher_dict = {
  'HDF5TPLT' : args.hdf5template,
  'BDACONF'  : args.bdaconfig,
  'NFILES'   : 10,
  'SYNCTIME' : time.time(),
  'INTTIME'  : 64 * 2048,
  'DISKMING' : 9999,
#  'SYNCTIME' : r['corr:feng_sync_time'],
#  'INTTIME'  : r['corr:acc_len'],
  'TRIGGER'  : 0,
}
  
# Reset various statistics counters

pubchan = 'hashpipe://%s/%d/set' % (args.host, 0)
for key, val in catcher_dict.iteritems():
   r.publish(pubchan, '%s=%s' % (key, val))
for v in ['NETWAT', 'NETREC', 'NETPRC']:
    r.publish(pubchan, '%sMN=99999' % (v))
    r.publish(pubchan, '%sMX=0' % (v))

print('Write baseline dist to redis')
# Use the bdaconfig file to write baselines per bin
baselines = {}
for n in range(4):
    baselines[n] = []

bdaconfig = np.loadtxt(config, dtype=np.int)
for i,t in enumerate(bdaconfig[:,2]):
    if (t==0): continue
    n = int(np.log(t)/np.log(2))
    if (n==4): n = 3
    baselines[n].append((bdaconfig[i,0], bdaconfig[i,1]))

for i in range(4):
    r.publish(pubchan, 'NBL%dSEC=%d'  % (2**(i+1), len(baselines[i])))

time.sleep(0.1)

# Release nethread hold
r.publish(pubchan, 'CNETHOLD=0')
r.publish(pubchan, 'TRIGGER=1')

