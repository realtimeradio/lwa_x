import time
import socket
import struct
import numpy

NANTS = 192
NCHANS = 8192 / 4 * 3
BYTES_PER_PACKET = 1024
MCNT_STEP = 2**18
NX = 8
TIME_DEMUX = 2
NT = 2

DESTIP = socket.gethostbyname("catcher")
DESTPORT = 10000

n_bls = (2*NANTS * (2*NANTS+1)) / 2 # 2xNANTS because that's how we deal with pols
n_bytes = NCHANS * n_bls * 2 * 4 # real/imag, bytes-per-word
n_pkts = (n_bytes / BYTES_PER_PACKET)
n_pkts_per_x = n_pkts / NX

print "MBytes per dump: %.2f" % (n_bytes / 1e6)
print "Packets per dump: %.4f" % (n_bytes / float(BYTES_PER_PACKET))
print "Packets per Xeng: %.4f" % n_pkts_per_x

s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

mcnt = 0
data = '\x00' * BYTES_PER_PACKET

while(True):
    start = time.time()
    offset = 0
    for pkt in range(n_pkts_per_x):
        #if pkt % 100000 == 0:
        #    print pkt
        for t in range(TIME_DEMUX):
            for xeng in range(NX):
                p = struct.pack('<QLHH', mcnt+NT*t, offset, xeng, BYTES_PER_PACKET) + data
                s.sendto(p, (DESTIP, DESTPORT))
        offset += BYTES_PER_PACKET
    elapsed = time.time() - start
    print "dump sent (MCNT: %d) in %.2fs (%.2f Gb/s)" % (mcnt, elapsed, TIME_DEMUX*n_bytes*8 / elapsed / 1e9)
    mcnt += MCNT_STEP

