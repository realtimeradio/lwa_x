#!/usr/bin/env python

import redis
import time
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Trigger data collection on the HERA catcher node',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('host', type=str, help='Host on which to capture data')
parser.add_argument('-r', dest='redishost', type=str, default='redishost', help='Host serving redis database')
parser.add_argument('-n', dest='nfiles', type=int, default=10, help='Number of files of data to capture')
parser.add_argument('-m', dest='msperfile', type=int, default=60000, help='Number of ms of data per file')
parser.add_argument('-t', dest='tag', type=str, default='none', help='A descriptive tag to go into data files')

args = parser.parse_args()

r = redis.Redis(args.redishost)

if len(args.tag) > 127:
  raise ValueError("Tag argument must be <127 characters!")

#Configure runtime parameters
catcher_dict = {
  'MSPERFIL' : args.msperfile,
  'NFILES'   : args.nfiles,
  'TAG'      : args.tag,
}
  
pubchan = 'hashpipe://%s/%d/set' % (args.host, 0)
for key, val in catcher_dict.iteritems():
   r.publish(pubchan, '%s=%s' % (key, val))

# Only trigger after the other parameters have had ample time to write
time.sleep(0.1)
r.publish(pubchan, "TRIGGER=1")
