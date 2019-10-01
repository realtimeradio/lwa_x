# BDA
# Check packet structure and contents in test mode

import numpy as np
import struct
import socket
import argparse 

Nbins = 5               # number of diff integration bins
Na = 192                # antennas
N_bl_per_block = 16384  # baselines within each block

parser = argparse.ArgumentParser(description='Test packet format and contents for BDA',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--host', type=str, default="127.0.0.1",
                     help='Host to capture data')
parser.add_argument('-p', dest='port', type=int, default=10000,
                     help='Host port to receive data')
parser.add_argument('-n', dest='n_inputs', type=int, default=704,
                     help='Nant*Npol that xGPU has been configured')
parser.add_argument('-v', dest='verbose', action='store_true', default=False,
                     help='Print some contents of each packet')
parser.add_argument('-c', dest='check', action='store_true', default=False,
                     help='Check the packet contents')

args = parser.parse_args()

# mcnt, bcnt, offset, ant0, ant1, xeng_id, payload_len
def decode_packet(packet):
   p = struct.unpack('>1Q2I4H1024i', packet)
   time, bcnt, offset, ant0, ant1, xeng_id, payload_len = p[0:7]
   data = np.asarray(p[7:])
   return time, bcnt, offset, ant0, ant1, xeng_id, payload_len, data

int_bin = {}
int_bin['baselines'] = {}
int_bin['data']      = {}

fakereal = 1
fakeimag = 2

for n in range(Nbins):
    int_bin['baselines'][n] = []

    # Ramp
    int_bin['data'][n] = (2**n)*np.repeat((np.arange(128)+fakereal),8)
    int_bin['data'][n][1::2] = -1*(2**n) * np.repeat((np.arange(128)+fakeimag),4)

    # Const
    #int_bin['data'][n] = (2**n)*np.ones(1024)
    #int_bin['data'][n][1::2] = -2*(2**n)

bdaconfig = np.loadtxt('../bda_config_192ants_nobda.txt', dtype=np.int)
for i,t in enumerate(bdaconfig[:,2]):
    if t==0: continue
    n = int(np.log(t)/np.log(2))
    int_bin['baselines'][n].append((bdaconfig[i,0], bdaconfig[i,1]))
#for a in range(Na):
#    int_bin['baselines'][0].append((a,a))


sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock.bind((args.host, args.port))

errors = 0
packets = 0

while True:
    try:
        pkt, addr = sock.recvfrom(4120) # buffer size is 1024 bytes
        packets += 1
        t,b,o,a0,a1,x,l,data = decode_packet(pkt)
        if args.verbose:
           print "{0:4d} {1:3d} {2:4d} {3:1d} {4:3d} {5:3d} {6:2d}".format(t, b//16384, b%16384, o, a0, a1, x), data[:8]
        if args.check:
           if (a0 > Na or a1 > Na): 
              print "Error! Received wrong antenna!"
           n = [y for y,v in int_bin['baselines'].items() if (a0,a1) in v][0]
           if not (np.all(data == int_bin['data'][n])):
              print "Error!", int_bin['data'][n][:32:8], o, n, data[:32:8]
              errors += 1
    except(KeyboardInterrupt):
       print("")
       print("Total number of packets captured: %d"%packets)
       print("Total number of errors: %d"%errors)
       break
