import socket
import numpy as np
import struct
import time
import cPickle as pickle


NANTS = 192
NCHANS = 8192 / 4 * 3 / 4
NXENG = 16
NWORDS = 2 * NANTS * (NANTS+1) / 2 * 4 * NCHANS

PAYLOAD_LEN = 4096
BYTES_PER_PACKET = PAYLOAD_LEN + 16
max_offset = NWORDS*8 / NXENG - PAYLOAD_LEN
print "Max offset:", max_offset

NPACKETS = NWORDS * 8 / PAYLOAD_LEN
print "Expecting %d packets" % NPACKETS

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
sock.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 80000000)
print sock.getsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF)

sock.bind(("10.80.40.251", 10000))

dout = np.zeros(2*NWORDS, dtype=np.int32)

n = 0
wait = True
start = time.time()
while True:
    data = sock.recv(BYTES_PER_PACKET)
    unpacked = struct.unpack('!QLHH', data[0:16])
    time_demux = (unpacked[0] % 4) >> 1
    timestamp = unpacked[0] - (unpacked[0] % 4)
    offset = unpacked[1] + (time_demux * NWORDS*4)
    xeng_id = unpacked[2]
    payload_len = unpacked[3]
    if wait: 
        start = time.time()
        this_window = timestamp
    if (not wait) or (offset == 0):
        n += 1
        wait = False
        #print timestamp, offset, xeng_id, payload_len
        #print dout[256*packet_num:256*(packet_num+1)].shape
        #print len(data[16:])
        #print dout.shape, offset/4, (offset+PAYLOAD_LEN)/4, dout[offset/4:(offset+PAYLOAD_LEN)/4].shape, np.fromstring(data[16:], dtype='>i').shape
        dout[offset/4:(offset+PAYLOAD_LEN)/4] = np.fromstring(data[16:], dtype='>i')
        if  (timestamp != this_window):
            print timestamp, this_window
            stop = time.time()
            print "%d packets received in %.2f seconds" % (n, stop - start)
            n = 0
            start = stop
            break


print "%d packets received in %.2f seconds" % (n, stop - start)

print "Dumping packet to disk"
with open("/tmp/packet.bin", "w") as fh:
    fh.write(dout.tostring())

