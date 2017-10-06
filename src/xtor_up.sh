#!/bin/bash

show_help_and_exit() {
cat <<.
Script to bring up the PAPER correlator.  Assumes that correlator is
completely down (e.g. after power up or after running xtor_down.sh).

Usage: xtor_up.sh [-f [-m MODE]] [-x|-c] [SLICES]

Passing the "-f" option brings up the F engines.
Passing the "-x" option brings up the X engines.
Passing the "-c" option brings up the CRC check engines.
Passing both -x and -c is not allowed.

Passing -m MODE allows setting the corner turner mode of the F engines.

MODE=0 == 8-way correlator (8 ROACH2s + 8 X Boxes) [default]
MODE=1 == 4-way correlator (4 ROACH2s + 4 X Boxes)
MODE=2 == 2-way correlator (2 ROACH2s + 2 X Boxes)
MODE=3 == 1-way correlator (1 ROACH2  + 1 X Boxes)

Examples:

  # Bring up all F and X engines of the PAPER correlator
  $ xtor_up.sh -f -x

  # Bring up slice 1 F engine only (pf1) with corner tuner mode 3
  $ xtor_up.sh -f -m 3 1

  # Bring up slices 1..4 CRC check engines only (px1..px4)
  $ xtor_up.sh -c 1..4
.
exit $1
}

do_f=
do_x=
do_c=
paper_init=paper_pktsock_init.sh
xc=X
ctmode=0

while getopts :fxcm:h-: opt
do
  case $opt in
    f)
      do_f=1
      ;;
    m)
      ctmode=$OPTARG
      ;;
    x)
      do_x=1
      ;;
    c)
      do_c=1
      paper_init=paper_crc_init.sh
      xc=C
      ;;
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

# If none of -f, -x, or -c were given, show help
if [ -z "${do_f}${do_x}${do_c}" ]
then
  show_help_and_exit
fi

# If both x and c were given, error out
if [ -n "${do_x}" -a -n "${do_c}" ]
then
  echo "error: cannot give both x and c options"
  exit 1
# Else, if neither f nor x nor c given, assume f and x
elif [ -z "${do_f}${do_x}${do_c}" ]
then
  do_f=1
  do_x=1
fi

# If no slices were given, do slices 1..8
if [ $# -eq 0 ]
then
  set {1..8}
fi

fhosts=$(echo "$@" | sed 's/^/pf/;s/ / pf/g')
xhosts=$(echo "$@" | sed 's/^/px/;s/ / px/g')

#echo do_f=${do_f}
#echo mode=${ctmode}
#echo do_x=${do_x}
#echo do_c=${do_c}
#echo "slices: ${@}"
#echo "fhosts: ${fhosts}"
#echo "xhosts: ${xhosts}"
#exit

if [ -n "$do_f" ]
then
  # Program ROACH2 FPGAs with roach2_fengine design
  echo Programming F Engines $fhosts
  for f in $fhosts; do adc16_init.rb -r 0x2a=0x6666 $f roach2_fengine; done

  # Initialize F engines
  echo Initializing F Engines $fhosts
  #paper_feng_init.rb -m $ctmode $fhosts
  hera_roach_feng_init.py $fhosts
fi

if [ -n "$do_x" -o "$do_c" ]
then
  # Start X engines
  echo Starting $xc Engine instances on $xhosts
  for x in $xhosts
  do
    ssh $x ${paper_init} 0 1 2 3 &
  done

  wait

  # Start hashpipe-redis gateways
  echo Starting hashpipe-redis gateways on $xhosts
  for x in $xhosts
  do
    # Start two gateways on each host: one gateway per NUMA node, two instances
    # per gateway.  The gateways run on the OS core (0 or 6)
    ssh $x "
      taskset 0x0001 hashpipe_redis_gateway.rb -g $x -i 0,1;
      taskset 0x0040 hashpipe_redis_gateway.rb -g $x -i 2,3
    "
  done

  # Let the gateways come up
  sleep 1

  # Turn off HOLD flag on all instances
  echo Enabling all X Engine network threads
  for i in 3 2 1 0
  do
    for x in $xhosts
    do
      redis-cli -h redishost publish hashpipe://$x/$i/set NETHOLD=0 > /dev/null
    done
  done

  sleep 1

  # Reset all NET{WAT,REC,PRC}M{N,X} counters
  echo Resetting all network stats counters
  for i in 3 2 1 0
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
fi
