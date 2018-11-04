ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS     =

AM_CPPFLAGS += -I"@HASHPIPE_INCDIR@"
AM_CPPFLAGS += -I"@XGPU_INCDIR@"
AM_CPPFLAGS += -I/usr/include/hdf5/serial/

# AM_CFLAGS is used for all C compiles
#AM_CFLAGS = -ggdb -fPIC -O3 -Wall -Werror -fno-strict-aliasing -mavx2
AM_CFLAGS = -fPIC -O3 -Wall -Werror -fno-strict-aliasing -mavx2 -march=native -funroll-loops

# Convenience variables to group source files
paper_databuf = paper_databuf.h \
	        paper_databuf.c

paper_threads = paper_fake_net_thread.c       \
		paper_gpu_cpu_output_thread.c \
		paper_gpu_output_thread.c     \
		paper_gpu_thread.c            \
		paper_fluff_thread.c          \
		hera_pktsock_thread.c         \
		paper_crc_thread.c            \
                hera_catcher_net_thread.c     \
                hera_catcher_disk_thread.c    \
                hera_gpu_output_thread.c      \
                #lzf/lzf/lzf_c.c \
                #lzf/lzf/lzf_d.c \
                #lzf/lzf_filter.c \
                #lzf_filter.h

paper_fluff = paper_fluff.h \
	      paper_fluff.c

hera_capture = hera_hdf5.h \
               hera_write_hdf5.c

# This is the paper_gpu plugin itself
lib_LTLIBRARIES = paper_gpu.la
paper_gpu_la_SOURCES  = $(paper_databuf) $(paper_threads) $(paper_fluff)
paper_gpu_la_LIBADD   = -lrt -lxgpu
paper_gpu_la_LDFLAGS  = -avoid-version -module -shared -export-dynamic
paper_gpu_la_LDFLAGS += -L"@XGPU_LIBDIR@" -Wl,-rpath,"@XGPU_LIBDIR@"
paper_gpu_la_LDFLAGS += -L"@HASHPIPE_LIBDIR@" -Wl,-rpath,"@HASHPIPE_LIBDIR@"
# Include HDF5 Libraries
paper_gpu_la_LDFLAGS += -L/usr/lib/x86_64-linux-gnu/hdf5/serial
paper_gpu_la_LDFLAGS += -lhdf5_hl -lhdf5
# Include HDF5 compression libraries
#paper_gpu_la_LDFLAGS += -L/usr/local/lib/python2.7/dist-packages/bitshuffle/plugin
#paper_gpu_la_LDFLAGS += -lh5LZF -lh5bshuf


# Installed scripts
dist_bin_SCRIPTS = ../scripts/paper_init.sh \
		   ../scripts/paper_crc_init.sh \
		   ../scripts/xtor_down.sh \
		   ../scripts/xtor_up.py \
                   ../scripts/catcher_init.sh \
                   ../scripts/catcher_up.py \
                   ../scripts/catcher_down.sh \
                   ../scripts/tweak-perf.sh \
                   ../scripts/hera_ctl.py \
                   ../scripts/hera_make_hdf5_template.py

# Non-installed test programs
noinst_PROGRAMS = paper_fluff_test hera_write_hdf5
paper_fluff_test_SOURCES = paper_fluff_test.c $(paper_fluff)
# The paper_fluff_test_CFLAGS is an autotools workaround.
# http://www.gnu.org/software/automake/manual/html_node/Objects-created-both-with-libtool-and-without.html
paper_fluff_test_CFLAGS = $(AM_CFLAGS)
paper_fluff_test_LDADD = -lrt

hera_write_hdf5_SOURCES = hera_write_hdf5 $(hera_capture)
# The hera_write_hdf5_CFLAGS is an autotools workaround.
# http://www.gnu.org/software/automake/manual/html_node/Objects-created-both-with-libtool-and-without.html
hera_write_hdf5_CFLAGS = -O3 -Wall -Werror -fno-strict-aliasing -mavx2 -march=native -funroll-loops
hera_write_hdf5_LDADD = -lrt -lpthread -lsz -lz -ldl -lm
hera_write_hdf5_LDFLAGS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a

# vi: set ts=8 noet :
