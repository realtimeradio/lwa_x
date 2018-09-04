#! /bin/bash

# Set high performance mode
for i in `ls /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor`; do echo performance > $i; done

# Set mtu
ifconfig eth3 mtu 9000

# Turn on pause requests
ethtool -A eth3 rx on

# Kernel buffer sizes
sysctl net.core.rmem_max=8388608
sysctl net.core.rmem_default=8388608

# Kill packets before the IP stack
iptables -t raw -A PREROUTING -i eth3 -p udp -j DROP

# Set interrupt coalescing
ethtool -C eth3 adaptive-rx off
ethtool -C eth3 rx-frames 8
ethtool -C eth3 rx-usecs 0

# Set ring sizes to max
ethtool -G eth3 rx 8192

# Set Receiver Side Steering
#ethtool -U eth3 flow-type udp4 src-ip 10.0.10.110 m 0.0.0.0 action 4 loc 1
#ethtool -U eth3 flow-type udp4 src-ip 10.0.10.111 m 0.0.0.0 action 5 loc 2
#ethtool -U eth3 flow-type udp4 src-ip 10.0.10.113 m 0.0.0.0 action 6 loc 3
#ethtool -U eth3 flow-type udp4 src-ip 10.0.10.114 m 0.0.0.0 action 7 loc 4
ethtool -U eth3 flow-type udp4 src-ip 10.80.40.1 m 255.255.255.255 action 4 loc 2


