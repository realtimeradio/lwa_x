#!/bin/bash

show_help_and_exit() {
cat <<.
Script to bring up the HERA correlator X-engines.  Assumes that correlator is
completely down (e.g. after power up or after running xtor_down.sh).

Usage: xtor_up.sh 
.
exit $1
}

paper_init=paper_init.sh
perf_tweaker=tweak-perf.sh
xc=X

while getopts :fxcm:h-: opt
do
  case $opt in
    h)
      show_help_and_exit 0
      ;;
    -)
      if [ $OPTARG == 'help' ]
      then
        show_help_and_exit 0
      else
        echo Invalid option: --$OPTARG
        show_help_and_exit 1
      fi
      ;;
    ?)
      echo Invalid option: -$OPTARG
      show_help_and_exit 1
      ;;
  esac
done
shift $((OPTIND-1))

# If no slices were given, do slices 1..16
if [ $# -eq 0 ]
then
  set {1..16}
fi

xhosts=$(echo "$@" | sed 's/^/px/;s/ / px/g')

# Run performance tweaking script
echo Running performance tweaking script on $xhosts
for x in $xhosts
do
  ssh root@$x ${perf_tweaker} &
done

wait

# Start X engines
echo Starting $xc Engine instances on $xhosts
for x in $xhosts
do
  ssh $x ${paper_init} 0 1 &
done

wait

# Start hashpipe-redis gateways
echo Starting hashpipe-redis gateways on $xhosts
for x in $xhosts
do
  # Start two gateways on each host: one gateway per NUMA node, one instance
  # per gateway.  The gateways run on the OS core (7 or 15)
  ssh $x "
    taskset 0x0080 hashpipe_redis_gateway.rb -g $x -i 0;
    taskset 0x8000 hashpipe_redis_gateway.rb -g $x -i 1;
  "
done

# Let the gateways come up
sleep 1

# Turn off HOLD flag on all instances
echo Enabling all X Engine network threads
for i in 1 0
do
  for x in $xhosts
  do
    redis-cli -h redishost publish hashpipe://$x/$i/set NETHOLD=0 > /dev/null
  done
done

sleep 1

# Reset all NET{WAT,REC,PRC}M{N,X} counters
echo Resetting all network stats counters
for i in 1 0
do
  for x in $xhosts
  do
    for k in NET{WAT,REC,PRC}
    do
      redis-cli -h redishost publish hashpipe://$x/$i/set ${k}MN=99999 > /dev/null
      redis-cli -h redishost publish hashpipe://$x/$i/set ${k}MX=0 > /dev/null
    done
  done
done
