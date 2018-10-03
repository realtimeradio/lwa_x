#! /bin/bash

# Set high performance mode
for i in `ls /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor`; do echo performance > $i; done

# Set mtu
ifconfig eth3 mtu 3200
ifconfig eth5 mtu 3200

# Turn on pause requests
ethtool -A eth3 rx on
ethtool -A eth5 rx on

# Kernel buffer sizes
sysctl net.core.rmem_max=8388608
sysctl net.core.rmem_default=8388608

# Kill packets before the IP stack
iptables -t raw -A PREROUTING -i eth3 -p udp -j DROP
iptables -t raw -A PREROUTING -i eth5 -p udp -j DROP

# Set interrupt coalescing
ethtool -C eth3 adaptive-rx off
ethtool -C eth3 rx-frames 8
ethtool -C eth3 rx-usecs 0
ethtool -C eth5 adaptive-rx off
ethtool -C eth5 rx-frames 8
ethtool -C eth5 rx-usecs 0

# Set ring sizes to max
ethtool -G eth3 rx 8192
ethtool -G eth5 rx 8192

# Set Receiver Side Steering
#ethtool -U eth3 flow-type udp4 src-ip 10.0.10.110 m 0.0.0.0 action 4 loc 1
#ethtool -U eth3 flow-type udp4 src-ip 10.0.10.111 m 0.0.0.0 action 5 loc 2
#ethtool -U eth3 flow-type udp4 src-ip 10.0.10.113 m 0.0.0.0 action 6 loc 3
#ethtool -U eth3 flow-type udp4 src-ip 10.0.10.114 m 0.0.0.0 action 7 loc 4
ethtool -U eth3 flow-type udp4 src-ip 10.0.10.110 m 0.0.0.0 action 4 loc 1
ethtool -U eth3 flow-type udp4 src-ip 10.0.10.111 m 0.0.0.0 action 4 loc 2
ethtool -U eth3 flow-type udp4 src-ip 10.0.10.113 m 0.0.0.0 action 5 loc 3
ethtool -U eth3 flow-type udp4 src-ip 10.0.10.114 m 0.0.0.0 action 5 loc 4
ethtool -U eth3 flow-type udp4 src-ip 10.0.10.201 m 0.0.0.0 action 6 loc 5
ethtool -U eth3 flow-type udp4 src-ip 10.0.10.109 m 0.0.0.0 action 6 loc 6
# Dump everything else on CPU7 (there shouldn't be much)
ethtool -U eth3 flow-type udp4 src-ip 10.0.10.1 m 255.255.255.255 action 4 loc 7

ethtool -U eth5 flow-type udp4 src-ip 10.0.10.110 m 0.0.0.0 action 12 loc 1
ethtool -U eth5 flow-type udp4 src-ip 10.0.10.111 m 0.0.0.0 action 12 loc 2
ethtool -U eth5 flow-type udp4 src-ip 10.0.10.113 m 0.0.0.0 action 13 loc 3
ethtool -U eth5 flow-type udp4 src-ip 10.0.10.114 m 0.0.0.0 action 13 loc 4
ethtool -U eth3 flow-type udp4 src-ip 10.0.10.201 m 0.0.0.0 action 14 loc 5
ethtool -U eth3 flow-type udp4 src-ip 10.0.10.109 m 0.0.0.0 action 14 loc 6
# Dump everything else on CPU7 (there shouldn't be much)
ethtool -U eth5 flow-type udp4 src-ip 10.0.10.1 m 255.255.255.255 action 15 loc 7


