import numpy as np
import struct
import socket
import time

Na = 352   # antennas
Nx = 16    # xeng
Nc = 328   # chan
Nt = 2     # demux
Ns = 4     # stokes
udp_ip = '10.20.1.66'
udp_port = 10000

bdaconfig = np.loadtxt('../bda_config.txt')

baselines = {}
for i in range(5):
    baselines[i] = []

for i,t in enumerate(bdaconfig[:,2]):
    baselines[int(np.log(t)/np.log(2))].append((bdaconfig[i,0],bdaconfig[i,1]))

for a0 in range(Na):
    baselines[1].append((a0,a0))

bcnt = 0
mcnt = 0
cntr = 0
offset = np.arange(3)
xeng_id = 0

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock.bind(('10.20.1.2',8500))
data = {}
for i in range(6):
    data[i] = (2**i)*np.ones(1024)
    data[i][1::2] = -1*(2**i)

# mcnt, bcnt, offset, ant0, ant1, xeng_id, payload_len
while True:
    cntr += 1
    mcnt = cntr*1024
    for a0,a1 in baselines[0]:
        for o in offset:
            pkt = struct.pack('>1Q2I4H1024i', mcnt, bcnt, o, a0, a1, xeng_id, 4096, *data[0])
            sock.sendto(pkt, (udp_ip, udp_port))
            time.sleep(0.001)
        bcnt += 1

    if (cntr%2 == 0):
        for a0,a1 in baselines[1]:
            for o in offset:
                pkt = struct.pack('>1Q2I4H1024i', mcnt, bcnt, o, a0, a1, xeng_id, 4096, *data[1])
                sock.sendto(pkt, (udp_ip, udp_port))
                time.sleep(0.001)
            bcnt += 1
      
    if (cntr%4 == 0):
        for a0,a1 in baselines[2]:
            for o in offset:
                pkt = struct.pack('>1Q2I4H1024i', mcnt, bcnt, o, a0, a1, xeng_id, 4096, *data[2])
                sock.sendto(pkt, (udp_ip, udp_port))
                time.sleep(0.001)
            bcnt += 1

    if (cntr%8 == 0):
        for a0,a1 in baselines[3]:
            for o in offset:
                pkt = struct.pack('>1Q2I4H1024i', mcnt, bcnt, o, a0, a1, xeng_id, 4096, *data[3])
                sock.sendto(pkt, (udp_ip, udp_port))
                time.sleep(0.001)
            bcnt += 1

    if (cntr%16 == 0):
        for a0,a1 in baselines[4]:
            for o in offset:
                pkt = struct.pack('>1Q2I4H1024i', mcnt, bcnt, o, a0, a1, xeng_id, 4096, *data[4])
                sock.sendto(pkt, (udp_ip, udp_port))
                time.sleep(0.001)
            bcnt += 1

