import numpy as np
import redis

r = redis.Redis('redishost')
pubchan = 'hashpipe://%s/%d/set' %('hera-sn1',0)

# Release nethread hold
r.publish(pubchan, 'CNETHOLD=0')
r.publish(pubchan,'TRIGGER=1')
