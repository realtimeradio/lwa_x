import time
import socket
import struct
import numpy

NANTS = 196
NCHANS = 8192 / 4 * 3
BYTES_PER_PACKET = 1024
MCNT_STEP = 2**18
NX = 8

DESTIP = "localhost"
DESTPORT = 10000

n_bls = (NANTS * (NANTS+1)) / 2
n_bytes = NCHANS * n_bls * 4 * 2 * 4 # stokes, real/imag, bytes-per-word
n_pkts = (n_bytes / BYTES_PER_PACKET)

print "MBytes per dump: %.2f" % (n_bytes / 1e6)
print "Packets per dump: %.4f" % (n_bytes / float(BYTES_PER_PACKET))

s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

mcnt = 0
data = '\x00' * BYTES_PER_PACKET

while(True):
    start = time.time()
    offset = 0
    for pkt in range(n_pkts):
        if pkt % 1000 == 0:
            #print pkt
            pass
        for xeng in range(NX):
            p = struct.pack('!QLHH', mcnt, offset, xeng, BYTES_PER_PACKET) + data
            #s.sendto(p, (DESTIP, DESTPORT))
        offset += BYTES_PER_PACKET
    elapsed = time.time() - start
    print "dump sent (MCNT: %d) in %.2fs (%.2f Gb/s)" % (mcnt, elapsed, n_bytes*8 / elapsed / 1e9)
    mcnt += MCNT_STEP

