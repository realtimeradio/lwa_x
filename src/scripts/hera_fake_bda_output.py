# Imitate 16 xengs, 2 times, all baselines 
# Send packets to catcher imitating the output
# of all the xengs in the chain.

import numpy as np
import argparse
import struct
import socket
import time

parser = argparse.ArgumentParser(description='Generate FAKE output to test catcher pipeline with baseline dependent averaging',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('config', type=str, help='BDA config file')
parser.add_argument('--catcher', type=str, default= '10.10.10.5',
                    help='IP address of the Catcher machine')
args = parser.parse_args()

Na = 352   # antennas
Nx = 16    # xeng
Nc = 384   # chan
Nt = 2     # demux
Ns = 4     # stokes
Nbins = 4  # number of diff integration bins
udp_ip = args.catcher #'10.80.40.251'
udp_port = 10000
inttimes = np.asarray([1,2,4,8,16])

int_bin = {}
int_bin['baselines'] = {}
int_bin['data']      = {}

for n in range(Nbins):
    int_bin['baselines'][n] = []
    int_bin['data'][n] = (2**n)*np.ones(1024, dtype=np.int32)
    int_bin['data'][n][1::2] = -2*(2**n)
    # convert to binary data for speedup
    int_bin['data'][n] = np.array(int_bin['data'][n]).byteswap().tostring()


# Populate baseline pairs from config file
bdaconfig = np.loadtxt(args.config, dtype=np.int)
for i,t in enumerate(bdaconfig[:,2]):
    if (t==0): continue
    n = int(np.log(t)/np.log(2))
    if (n==4): n = 3
    int_bin['baselines'][n].append((bdaconfig[i,0], bdaconfig[i,1]))

for b in int_bin['baselines'].keys():
    print b,len(int_bin['baselines'][b])

bcnt = 0; mcnt = 0; ctr = 0;
sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
while True:
   ctr += 2
   mcnt = int(500e6 * ctr / (2 * 8192))
   for ns in np.logspace(1, Nbins, num=Nbins, base=2, dtype=np.int):
       if (ctr%ns == 0):
           nb = int(np.log2(ns)) - 1
           print 'Sending: %d Baselines: %d Bcnt: %d' % (ns, len(int_bin['baselines'][nb]), bcnt)
           for (a0,a1) in int_bin['baselines'][nb]:
               for xeng_id in range(Nx*Nt):
                   b = mcnt+((xeng_id//Nx)*2)
                   pkt = struct.pack('>1Q2I4H', b, bcnt, 0, a0, a1, xeng_id, 4096) + int_bin['data'][nb]
                   sock.sendto(pkt, (udp_ip, udp_port))
                   time.sleep(1e-3)
                   pkt = struct.pack('>1Q2I4H', b, bcnt, 1, a0, a1, xeng_id, 4096) + int_bin['data'][nb]
                   sock.sendto(pkt, (udp_ip, udp_port))
                   time.sleep(1e-3)
                   pkt = struct.pack('>1Q2I4H', b, bcnt, 2, a0, a1, xeng_id, 4096) + int_bin['data'][nb]
                   sock.sendto(pkt, (udp_ip, udp_port))
                   time.sleep(1e-3)
               bcnt += 1

#def send_xeng(xeng_id):
#    print 'Xeng: ', xeng_id
#    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
#    #sock.bind(("10.10.10.5",10000+xeng_id))
#    bcnt = 0; mcnt = 0; ctr = 1;
#    while True:
#        for nb in np.nonzero(np.logical_not(ctr%inttimes))[0]:
#            for (a0,a1) in int_bin['baselines'][nb]:
#                b = mcnt+((xeng_id//Nx)*2)
#                pkt = struct.pack('>1Q2I4H', b, bcnt, 0, a0, a1, xeng_id, 4096) + int_bin['data'][nb]
#                sock.sendto(pkt, (udp_ip, udp_port))
#                pkt = struct.pack('>1Q2I4H', b, bcnt, 1, a0, a1, xeng_id, 4096) + int_bin['data'][nb]
#                sock.sendto(pkt, (udp_ip, udp_port))
#                pkt = struct.pack('>1Q2I4H', b, bcnt, 2, a0, a1, xeng_id, 4096) + int_bin['data'][nb]
#                sock.sendto(pkt, (udp_ip, udp_port))
#            bcnt += 1
#        ctr  += 1
#

#threads = []
#for xeng_id in range(16): 
#    t = threading.Thread(target=send_xeng, args=(xeng_id,))
#    threads.append(t)
#    t.start()

#xeng_id = 0
#send_xeng(xeng_id)

