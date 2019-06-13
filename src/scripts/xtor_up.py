#!/usr/bin/env python

import redis
import time
import argparse
import subprocess

perf_tweaker = 'tweak-perf.sh'
paper_init = 'paper_init.sh'
paper_init_ibv = 'paper_init_ibv.sh'

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

parser.add_argument('hosts', type=str, nargs='+', help='Hosts to intialize')
parser.add_argument('-r', dest='redishost', type=str, default='redishost', help='Host serving redis database')
parser.add_argument('-t', dest='timeslices', type=int, default=2,
                    help='Number of independent correlators. E.g. 2 => Even/odd correlator')
parser.add_argument('-i', dest='ninstances', type=int, default=2,
                    help='Number of pipeline instances per host')
parser.add_argument('--runtweak', dest='runtweak', action='store_true', default=False,
                    help='Run the tweaking script %s on X-hosts prior to starting the correlator' % perf_tweaker)
parser.add_argument('--ibverbs', dest='ibverbs', action='store_true', default=False,
                    help='Use the IB Verbs netthread. Experimental!')

args = parser.parse_args()
hosts = args.hosts # Too lazy to keey typing this
nhosts = len(hosts)
nhosts_per_timeslice = nhosts / args.timeslices

assert args.ninstances == 2, 'Sorry, anything other than ninstances=2 is not supported!'

r = redis.Redis(args.redishost)

# Run performance tweaking script
if args.runtweak:
    run_on_hosts(hosts, perf_tweaker, user='root', wait=True)

# Start X-Engines
if args.ibverbs:
    run_on_hosts(hosts, [paper_init_ibv, '0', '1'], wait=True) # two instances per host
else:
    run_on_hosts(hosts, [paper_init, '0', '1'], wait=True) # two instances per host

# Start hashpipe<->redis gateways
cpu_masks = ['0x0080', '0x8000']
for host in hosts:
    for i in range(args.ninstances):
        run_on_hosts([host], ['taskset', cpu_masks[i], 'hashpipe_redis_gateway.rb', '-g', host, '-i', '%d'%i])

# Wait for the gateways to come up
time.sleep(3)

# Configure the X-engines as even/odd correlators
for hn, host in enumerate(hosts):
   for i in range(args.ninstances):
      key = 'hashpipe://%s/%d/set' % (host, i)
      r.publish(key, 'TIMEIDX=%d' % (hn//nhosts_per_timeslice))

# Let the network threads begin processing
for hn, host in enumerate(hosts):
   for i in range(args.ninstances):
      key = 'hashpipe://%s/%d/set' % (host, i)
      r.publish(key, 'NETHOLD=0')

time.sleep(2)

# Reset various statistics counters
for hn, host in enumerate(hosts):
   for i in range(args.ninstances):
      key = 'hashpipe://%s/%d/set' % (host, i)
      for v in ['NETWAT', 'NETREC', 'NETPRC']:
          r.publish(key, '%sMN=99999' % (v))
          r.publish(key, '%sMX=0' % (v))
