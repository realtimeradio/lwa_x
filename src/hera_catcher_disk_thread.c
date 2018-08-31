/*
 * hera_catcher_disk_thread.c
 *
 * Writes correlator to disk as hdf5 files
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/sendfile.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <hdf5.h>

#include "hashpipe.h"
#include "paper_databuf.h"

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

#define N_DATA_DIMS (4)
/*
Create an extensible dataset for visdata, visdata_diff, flags, nsamples
see https://gist.github.com/simleb/5205083/
*/
void make_extensible_hdf5(hid_t id)
{
    hsize_t dims[N_DATA_DIMS] = {0, N_CHAN_TOTAL, VIS_MATRIX_ENTRIES_PER_CHAN, N_STOKES};
    hsize_t max_dims[N_DATA_DIMS] = {H5S_UNLIMITED, N_CHAN_TOTAL, VIS_MATRIX_ENTRIES_PER_CHAN, N_STOKES};
    
    hid_t file_space = H5Screate_simple(N_DATA_DIMS, dims, max_dims);
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    hsize_t chunk_dims[N_DATA_DIMS] = {2, N_CHAN_TOTAL, VIS_MATRIX_ENTRIES_PER_CHAN, N_STOKES};
    H5Pset_chunk(plist, N_DATA_DIMS, chunk_dims);
    
    // Now we have the dataspace properties, create the datasets
    H5Dcreate(id, "visdata",      H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Dcreate(id, "visdata_diff", H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Dcreate(id, "nsamples",     H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Dcreate(id, "flags",        H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(file_space);
}
/*
Turn an mcnt into a UNIX time in double-precision.
*/
double mcnt2time(uint32_t mcnt, uint32_t sync_time)
{
    return (sync_time + mcnt) * (N_CHAN_TOTAL_GENERATED / (double)FENG_SAMPLE_RATE);
}

hid_t open_hdf5_from_template(char * sourcename, char * destname)
{
    int read_fd, write_fd;
    struct stat stat_buf;
    off_t offset = 0;
    
    read_fd = open(sourcename, O_RDONLY);
    if (read_fd < 0) {
        return -1;
    }

    // Get the file size
    fstat(read_fd, &stat_buf);

    write_fd = open(destname, O_WRONLY | O_CREAT, stat_buf.st_mode);
    if (write_fd < 0) {
        return -1;
    }

    sendfile(write_fd, read_fd, &offset, stat_buf.st_size);
    
    close(read_fd);
    close(write_fd);
    return H5Fopen(sourcename, H5F_ACC_RDWR, H5P_DEFAULT);
}

static void *run(hashpipe_thread_args_t * args)
{
    // Local aliases to shorten access to args fields
    // Our input buffer is a paper_input_databuf
    // Our output buffer is a paper_gpu_input_databuf
    hera_catcher_input_databuf_t *db_in = (hera_catcher_input_databuf_t *)args->ibuf;
    hashpipe_status_t st = args->st;
    const char * status_key = args->thread_desc->skey;
     
    // Buffers for file name strings
    char template_fname[128];
    char hdf5_fname[128];

    // Variables for sync time and computed gps time / JD
    uint32_t sync_time = 0;
    uint32_t gps_time;
    uint32_t julian_time;

    // How many integrations to dump to a file before starting the next one
    uint32_t ms_per_file;

    // Get template filename from redis
    hashpipe_status_lock_safe(&st);
    hgets(st.buf, "HDF5TPLT", 128, template_fname);
    
    // Get time that F-engines were last sync'd
    hgetu4(st.buf, "SYNCTIME", &sync_time);

    hgetu4(st.buf, "MSPERFILE", &ms_per_file);

    // Init status variables
    hputi8(st.buf, "DISKMCNT", 0);
    hashpipe_status_unlock_safe(&st);

    /* Loop */
    int32_t *db_in32;
    int rv;
    int curblock_in=0;
    int chan_index, bl;
    int32_t *data_even;
    int32_t *data_odd;
    double curr_file_time = -1.0;
    hid_t hdf5_id, header_id, extra_header_id, data_group_id, dataset_id;
    hid_t visdata_id, visdata_diff_id, flags_id, nsamples_id, mem_space;
    double file_start_t, file_stop_t, file_duration;
    int64_t file_obs_id, file_nblts=0;
    float gbps, min_gbps;

    struct timespec start, finish;
    
    // Allocate the buffers to be used for summing channels
    data_even  = (int32_t *)malloc(VIS_MATRIX_ENTRIES_PER_CHAN * 2);
    data_odd   = (int32_t *)malloc(VIS_MATRIX_ENTRIES_PER_CHAN * 2);
    // Define the memory space used by these buffers for HDF5 access
    hsize_t dims[N_DATA_DIMS] = {1, 1, VIS_MATRIX_ENTRIES_PER_CHAN, N_STOKES};
    mem_space = H5Screate_simple(N_DATA_DIMS, dims, NULL);
    mem_space = mem_space; //FIXME

    while (run_threads()) {
        // Note waiting status,
        // query integrating status
        // and, if armed, start count
        hashpipe_status_lock_safe(&st);
        hputs(st.buf, status_key, "waiting");
        hashpipe_status_unlock_safe(&st);

        // Wait for new input block to be filled
        while ((rv=hera_catcher_input_databuf_wait_filled(db_in, curblock_in)) != HASHPIPE_OK) {
            if (rv==HASHPIPE_TIMEOUT) {
                hashpipe_status_lock_safe(&st);
                hputs(st.buf, status_key, "blocked_in");
                hashpipe_status_unlock_safe(&st);
                continue;
            } else {
                hashpipe_error(__FUNCTION__, "error waiting for filled databuf");
                pthread_exit(NULL);
                break;
            }
        }

        // Got a new data block, update status
        hashpipe_status_lock_safe(&st);
        hputs(st.buf, status_key, "writing");
        hputi4(st.buf, "DISKBKIN", curblock_in);
        hputu8(st.buf, "DISKMCNT", db_in->block[curblock_in].header.mcnt);
        hashpipe_status_unlock_safe(&st);
        
        db_in32 = (int32_t *)db_in->block[curblock_in].data;

        clock_gettime(CLOCK_MONOTONIC, &start);
        gps_time = mcnt2time(db_in->block[curblock_in].header.mcnt, sync_time);
        julian_time = 2440587.5 + (gps_time / (double)(86400.0));
        if ((gps_time - curr_file_time) > ms_per_file) {
            // If a file is open:
            if (curr_file_time < 0) {
                // Write the meta-data we now know at file close
                dataset_id = H5Dopen(extra_header_id, "stopt", H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_stop_t);
                H5Dclose(dataset_id);
                dataset_id = H5Dopen(extra_header_id, "duration", H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_duration);
                H5Dclose(dataset_id);
                dataset_id = H5Dopen(header_id, "Nblts", H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_nblts);
                // Close dataset
                H5Dclose(dataset_id);
                // Close groups
                H5Gclose(header_id);
                H5Gclose(extra_header_id);
                H5Gclose(data_group_id);
                // Close file
                H5Fclose(hdf5_id);
            }
            // And now start a new file
            curr_file_time = gps_time;
            file_nblts = 0;
            file_start_t = gps_time;
            file_obs_id = (int64_t)gps_time;
            sprintf(hdf5_fname, "zen.%7.5d.uvh5", julian_time);
            hdf5_id = open_hdf5_from_template(hdf5_fname, template_fname);
            // Open HDF5 groups
            header_id = H5Gopen(hdf5_id, "Header", H5P_DEFAULT);
            extra_header_id = H5Gopen(header_id, "extra_keywords", H5P_DEFAULT);
            data_group_id = H5Gopen(hdf5_id, "Data", H5P_DEFAULT);
            make_extensible_hdf5(data_group_id);
            // Get IDs for all the extensible datasets
            visdata_id      = H5Dopen(data_group_id, "visdata", H5P_DEFAULT);
            visdata_diff_id = H5Dopen(data_group_id, "visdata_diff", H5P_DEFAULT);
            flags_id        = H5Dopen(data_group_id, "flags", H5P_DEFAULT);
            nsamples_id     = H5Dopen(data_group_id, "nsamples", H5P_DEFAULT);
visdata_id      = visdata_id;
visdata_diff_id = visdata_diff_id;
flags_id        = flags_id       ;          
nsamples_id     = nsamples_id    ;
            // Write meta-data values we know at file-open
            dataset_id = H5Dopen(extra_header_id, "obs_id", H5P_DEFAULT);
            H5Dwrite(dataset_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_obs_id);
            H5Dclose(dataset_id);
            dataset_id = H5Dopen(extra_header_id, "startt", H5P_DEFAULT);
            H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_start_t);
            H5Dclose(dataset_id);
            
            
        }
        file_stop_t = gps_time;
        file_duration = file_stop_t - file_start_t;
        file_nblts += VIS_MATRIX_ENTRIES_PER_CHAN;
            
        // Write lst_array, time_array, uvw_array, zenith_dec, zenith_ra
        // Write data!
        // First accumulate multiple channels. Then make the even/odd sum / differences
        
        for(chan_index=0; chan_index<N_CHAN_TOTAL; chan_index+=1){
            if ((chan_index % CATCHER_CHAN_SUM) == 0) {
                // Start new accumulation (count to VIS_MATRIX_ENTRIES_PER_CHAN*2 for real/imag)
                for(bl=0; bl<VIS_MATRIX_ENTRIES_PER_CHAN*2; bl+=1){
                    data_even[bl] = db_in32[bl + 2*(VIS_MATRIX_ENTRIES_PER_CHAN*chan_index)];
                    data_odd[bl]  = db_in32[bl + 2*(VIS_MATRIX_ENTRIES_PER_CHAN*chan_index + VIS_MATRIX_ENTRIES)];
                }
            } else {
                // Add to existing accumulation
                for(bl=0; bl<VIS_MATRIX_ENTRIES_PER_CHAN*2; bl+=1){
                    data_even[bl] += db_in32[bl + 2*(VIS_MATRIX_ENTRIES_PER_CHAN*chan_index)];
                    data_odd[bl]  += db_in32[bl + 2*(VIS_MATRIX_ENTRIES_PER_CHAN*chan_index + VIS_MATRIX_ENTRIES)];
                }
            }
            if ((chan_index % CATCHER_CHAN_SUM) == (CATCHER_CHAN_SUM-1)) {
                // Write to file
                for(bl=0; bl<VIS_MATRIX_ENTRIES_PER_CHAN*2; bl+=1){
                data_even[bl] = data_even[bl];
                data_odd[bl] = data_odd[bl];
                }
            }
        }
        
        
            

        clock_gettime(CLOCK_MONOTONIC, &finish);

        // Note processing time
        hashpipe_status_lock_safe(&st);
        // Bits per fluff / ns per fluff = Gbps
        hgetr4(st.buf, "DISKMING", &min_gbps);
        gbps = (float)(8L*VIS_MATRIX_ENTRIES*2*4)/ELAPSED_NS(start,finish); //Gigabits / s
        hputr4(st.buf, "DISKGBPS", gbps);
        if(min_gbps == 0 || gbps < min_gbps) {
          hputr4(st.buf, "DISKMING", gbps);
        }
        hashpipe_status_unlock_safe(&st);

        // Mark input block as free and advance
        hera_catcher_input_databuf_set_free(db_in, curblock_in);
        curblock_in = (curblock_in + 1) % db_in->header.n_block;

        /* Check for cancel */
        pthread_testcancel();
    }

    // Thread success!
    return NULL;
}

static hashpipe_thread_desc_t disk_thread = {
    name: "hera_catcher_disk_thread",
    skey: "DISKSTAT",
    init: NULL,
    run:  run,
    ibuf_desc: {hera_catcher_input_databuf_create},
    obuf_desc: {NULL}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&disk_thread);
}
