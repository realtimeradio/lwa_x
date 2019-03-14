# BDA
# Check packet structure and contents in test mode

import numpy as np
import struct
import socket
import argparse 

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

def tri_index(i,j):
    return (i* (i+1))/2 + j

def baseline_index(n):
   """ Generate a baseline pair to index map given
       number of inputs 'n' """
   nant_2 = (n/2) // 2
   triangle_size = ((nant_2 + 1) * nant_2)/2
   middle_rect_offset = triangle_size
   last_cell_offset = 4 * middle_rect_offset - nant_2 - 1
   
   idx_map = {}
   
   for a0 in range(0, n/2, 1):
      for a1 in range(a0, n/2, 1):
          delta = a1-a0
          if (delta > nant_2):
              cell_index = last_cell_offset - tri_index(nant_2-2-a0, (n/2)-1-a1)
          elif (a1 < nant_2):
              cell_index = tri_index(a1,a0)
          else:
              cell_index = middle_rect_offset + (a1-nant_2)*(nant_2+1) + (nant_2-delta) 
          idx_map[cell_index] = (a0,a1)
   return idx_map

# mcnt, bcnt, offset, ant0, ant1, xeng_id, payload_len
def decode_packet(packet):
   p = struct.unpack('>1Q2I2H1025i', packet)
   time, bcnt, offset, ant0, ant1, xeng_id, payload_len = p[0:7]
   data = np.asarray(p[7:])
   return time, bcnt, offset, ant0, ant1, xeng_id, payload_len, data


idx_map = baseline_index(args.n_inputs)

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
           print t, b//4096, b%4096, o, a0, a1, data[:8], data[-8:]
        if args.check:
           if not (np.all(data[::2] == 1)):
              print "Error! Real", idx_map[b], o, data[::2]
              errors += 1
           if not (np.all(data[1::2] == -1)):
              print "Error! Imag", idx_map[b], o, data[1::2]
              errors += 1
    except(KeyboardInterrupt):
       print("")
       print("Total number of packets captured: %d"%packets)
       print("Total number of errors: %d"%errors)
       break
