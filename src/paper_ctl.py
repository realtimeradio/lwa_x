#!/usr/bin/env python

import redis
import time
import argparse

def mcnts_per_second(sample_rate=500e6, spectra_len=8192, round_to=2):
    """
    Return number of MCNTs in 1 second. For HERA, but not in general,
    this is the same as the number of spectra in 1 second.
    sample_rate: ADC clock rate in MHz
    spectra_len: Number of frequency channels in 1 F-Engine spectra (prior to any subselecting)
    round_to   : Output the number of mcnts per second rounded up to a multiple of this number
    """
    return round_to * (int(sample_rate / (spectra_len * 2) + (round_to - 1)) // round_to)


parser = argparse.ArgumentParser(description='Turn on and off the HERA correlator',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('action',type=str,
                    help = 'Action: "start"|"stop": start|stop the correlator.')
parser.add_argument('-n', dest='acclen', type=int, default=64*2048,
                    help ='Number of spectra to accumulate')
parser.add_argument('-t', dest='starttime', type=float, default=time.time()+5,
                    help ='UNIX time to start observations. Default is NOW+5s')
parser.add_argument('-x', dest='xhosts', type=int, default=16,
                    help ='Number of GPU hosts')
parser.add_argument('-i', dest='xinstances_per_host', type=int, default=2,
                    help ='Number of X-engine instances per GPU host')
parser.add_argument('-r', dest='redishost', type=str, default='redishost',
                    help ='Hostname of redis server')
args = parser.parse_args()

if args.action not in ['start', 'stop']:
    print 'Available actions are "start" and "stop"'
    exit

# Make status keys for the various X-engine instances which identify their redis entries
status_keys = []
xinstances = args.xhosts * args.xinstances_per_host
for xhost in range(1, args.xhosts+1):
    for instance in range(args.xinstances_per_host):
        status_keys += ['hashpipe://px%d/%d/status' % (xhost, instance)]

rdb = redis.Redis(args.redishost)

if args.action == 'start':
    gpumcnts = [rdb.hget(skey, 'GPUMCNT') for skey in status_keys]
    
    n_missing_mcnts = gpumcnts.count(None)
    if n_missing_mcnts == xinstances:
        print 'Could\'t find any GPUMCNT values. Cannot start'
        exit
    elif n_missing_mcnts > 0:
        print 'WARNING: Missing GPUMCNT values from instances' % n_missing_mcnts

    # try and convert these to integers. Skip on failure (probably means the value is "None")
    for i in range(len(gpumcnts)):
        try:
            gpumcnts[i] = int(gpumcnts[i])
        except:
            pass

    # Delay before starting, in MCNTs, rounded to 4096
    time_delay = args.starttime - time.time()
    delay_mcnts = 4096 * ((int(time_delay * mcnts_per_second()) + 2047) / 4096)
    
    # MCNT to start on. Don't need to round, becase gpumcnts has to be aligned to a block boundary
    trig_mcnt = max(gpumcnts) + delay_mcnts

    print 'Current MCNTs are:', gpumcnts
    print 'MCNTs per second: %.1f' % mcnts_per_second()
    print 'Trigger MCNT: %d' % trig_mcnt
    print 'Trigger time is %.1f seconds in the future' % time_delay

    # Use the hashpipe publish channel to update keys in all status buffers.
    # See the docstring at https://github.com/david-macmahon/rb-hashpipe/blob/master/bin/hashpipe_redis_gateway.rb
    # for details about formatting
    
    msg = 'INTSYNC=%d\nINTCOUNT=%d\nINTSTAT=start\nOUTDUMPS=0' % (trig_mcnt, args.acclen)
    rdb.publish("hashpipe:///set", msg)

