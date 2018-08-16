#!/bin/bash

show_help_and_exit() {
cat <<.
Script to bring up the HERA correlator data acquisition server.
.
exit $1
}

init=catcher_init.sh
catcherhost=hera-sn1

# Start X engines
echo Starting catcher instance on $catcherhost
echo "ssh $catcherhost \"${init} 0\" &"
ssh $catcherhost "${init} 0" &
wait

# Start hashpipe-redis gateways
echo Starting hashpipe-redis gateways on $catcherhost
echo "ssh $catcherhost \"taskset 0x0001 hashpipe_redis_gateway.rb -g $catcherhost -i 0\""
ssh $catcherhost "taskset 0x0001 hashpipe_redis_gateway.rb -g $catcherhost -i 0"

# Let the gateways come up
sleep 1

# Turn off HOLD flag on all instances
echo Enabling catcher network threads
redis-cli -h redishost publish hashpipe://$catcherhost/0/set CNETHOLD=0 > /dev/null

sleep 1

# Reset all NET{WAT,REC,PRC}M{N,X} counters
echo Resetting all network stats counters
for k in NET{WAT,REC,PRC}
do
  redis-cli -h redishost publish hashpipe://$catcherhost/0/set ${k}MN=99999 > /dev/null
  redis-cli -h redishost publish hashpipe://$catcherhost/0/set ${k}MX=0 > /dev/null
done
