#!/usr/bin/env python

import os
import redis
import time
import argparse
import subprocess

perf_tweaker = 'tweak-perf-sn.sh'
init = 'hera_catcher_init.sh'
python_source_cmd = ['source', '~/hera-venv/bin/activate']
template_cmd = ['hera_make_hdf5_template.py']

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

parser = argparse.ArgumentParser(description='Start the HERA X-engines',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('host', type=str, help='Host to intialize')
parser.add_argument('-r', dest='redishost', type=str, default='redishost', help='Host serving redis database')
parser.add_argument('-t', dest='hdf5template', type=str, default='/tmp/template.h5', 
                    help='Place to put HDF5 header template file')
parser.add_argument('--runtweak', dest='runtweak', action='store_true', default=False,
                    help='Run the tweaking script %s on X-hosts prior to starting the correlator' % perf_tweaker)
parser.add_argument('--redislog', dest='redislog', action='store_true', default=False,
                    help='Use the redis logger to duplicate log messages on redishost\'s log-channel pubsub stream')
parser.add_argument('--pypath', dest='pypath', type=str, default="/home/hera/hera-venv",
                    help='The path to a python virtual environment which will be activated prior to running paper_init. ' +
                         'Only relevant if using the --redislog flag, which uses a python redis interface')

args = parser.parse_args()

r = redis.Redis(args.redishost)

# Run performance tweaking script
if args.runtweak:
    run_on_hosts([args.host], perf_tweaker, user='root', wait=True)

# Start Catcher
if args.redislog:
    python_source_cmd = ["source", os.path.join(args.pypath, "bin/activate")+";"]
    run_on_hosts([args.host], python_source_cmd + ['cd', '/data;', init, '-r', '0'], wait=True)
else:
    run_on_hosts([args.host], ['cd', '/data;', init, '0'], wait=True)
time.sleep(15)

# Start hashpipe<->redis gateways
cpu_mask = '0x0004'
run_on_hosts([args.host], ['taskset', cpu_mask, 'hashpipe_redis_gateway.rb', '-g', args.host, '-i', '0'])

# Wait for the gateways to come up
time.sleep(15)

# Generate the meta-data template
run_on_hosts([args.host], python_source_cmd + [';'] + template_cmd + ['-c', '-r', args.hdf5template], wait=True)

#Configure runtime parameters
catcher_dict = {
  'HDF5TPLT' : args.hdf5template,
  'NFILES'   : 0,
  'SYNCTIME' : r['corr:feng_sync_time'],
  'INTTIME'  : r['corr:acc_len'],
  'TRIGGER'  : 0,
}
  
# Reset various statistics counters

pubchan = 'hashpipe://%s/%d/set' % (args.host, 0)
for key, val in catcher_dict.iteritems():
   r.publish(pubchan, '%s=%s' % (key, val))
for v in ['NETWAT', 'NETREC', 'NETPRC']:
    r.publish(pubchan, '%sMN=99999' % (v))
    r.publish(pubchan, '%sMX=0' % (v))

# Release nethread hold
r.publish(pubchan, 'CNETHOLD=0')
