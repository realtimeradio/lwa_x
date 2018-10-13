#!/bin/bash

# Add directory containing this script to PATH
PATH="$(dirname $0):${PATH}"

hostname=`hostname -s`

function getip() {
  out=$(host $1) && echo $out | awk '{print $NF}'
}

myip=$(getip $(hostname))

function init() {
  instance=0
  mask=0x00ff
  bindhost=eth3
  netcpu=0
  outcpu=1

  echo taskset $mask \
  hashpipe -p paper_gpu -I $instance \
    -o BINDHOST=$bindhost \
    -c $netcpu hera_catcher_net_thread \
    -c $outcpu hera_catcher_disk_thread

  taskset $mask \
  hashpipe -p paper_gpu -I $instance \
    -o BINDHOST=$bindhost \
    -c $netcpu hera_catcher_net_thread \
    -c $outcpu hera_catcher_disk_thread \
     < /dev/null \
    1> ~/catcher.out.$instance \
    2> ~/catcher.err.$instance &
}

if [ -z "$1" ]
then
  echo "Usage: $(basename $0) INSTANCE_ID [...]"
  exit 1
fi

for instidx in "$@"
do
  echo Starting instance catcher/$instidx
  init
  echo Instance catcher/$instidx pid $!
  # Sleep to let instance come up
  sleep 10
done

# Zero out MISSEDPK counts
for instidx in "$@"
do
  echo Resetting MISSEDPK counts for catcher/$instidx
  hashpipe_check_status -I $instidx -k MISSEDPK -s 0
done
