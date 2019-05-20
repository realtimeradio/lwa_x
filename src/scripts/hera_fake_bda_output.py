# Imitate 16 xengs, 2 times, all baselines 
# Send packets to catcher imitating the output
# of all the xengs in the chain.

import numpy as np
import struct
import socket
import time

Na = 352   # antennas
Nx = 16    # xeng
Nc = 328   # chan
Nt = 2     # demux
Ns = 4     # stokes
Nbins = 5  # number of diff integration bins
udp_ip = '10.80.40.251'
udp_port = 10000

int_bin = {}
int_bin['baselines'] = {}
int_bin['data']      = {}

for n in range(Nbins):
    int_bin['baselines'][n] = []
    int_bin['data'][n] = (2**n)*np.ones(1024)
    int_bin['data'][n][1::2] = -1*(2**n)

# Populate baseline pairs from config file
bdaconfig = np.loadtxt('../bda_config.txt', dtype=np.int)
for i,t in enumerate(bdaconfig[:,2]):
    n = int(np.log(t)/np.log(2))
    int_bin['baselines'][n].append((bdaconfig[i,0], bdaconfig[i,1]))
for a in range(Na):
    int_bin['baselines'][0].append((a,a))

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

bcnt = 0; mcnt = 0; ctr = 0
offset = np.arange(3)

while True:
    ctr += 2
    mcnt = ctr 
    for ns in np.logspace(0, Nbins-1, num=Nbins, base=2, dtype=np.int):
        if (ctr%ns == 0):
            print 'Sending: %d'%ns
            nb = int(np.log(ns)/np.log(2))
            for (a0,a1) in int_bin['baselines'][nb]:
                for xeng_id in range(Nx*Nt):
                    pkt = struct.pack('>1Q2I4H1024i', mcnt+((xeng_id//Nx)*2), bcnt, 0, a0, a1, xeng_id, 4096, *int_bin['data'][nb]);
                    sock.sendto(pkt, (udp_ip, udp_port))
                    time.sleep(0.000001)
                    pkt = struct.pack('>1Q2I4H1024i', mcnt+((xeng_id//Nx)*2), bcnt, 1, a0, a1, xeng_id, 4096, *int_bin['data'][nb]);
                    sock.sendto(pkt, (udp_ip, udp_port))
                    time.sleep(0.000001)
                    pkt = struct.pack('>1Q2I4H1024i', mcnt+((xeng_id//Nx)*2), bcnt, 2, a0, a1, xeng_id, 4096, *int_bin['data'][nb]);
                    sock.sendto(pkt, (udp_ip, udp_port))
                    time.sleep(0.000001)
                bcnt += 1 
