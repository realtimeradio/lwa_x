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
#include <smmintrin.h>
#include <immintrin.h>

#include "hashpipe.h"
#include "paper_databuf.h"
#include "lzf_filter.h"

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

#define N_DATA_DIMS (4)
#define N_CHAN_PROCESSED (N_CHAN_TOTAL / CATCHER_CHAN_SUM)
#define N_BL_PER_WRITE (8)

static hid_t complex_id;

typedef struct {
    hid_t file_id;
    hid_t header_gid;
    hid_t data_gid;
    hid_t extra_keywords_gid;
    hid_t visdata_did;
    hid_t flags_did;
    hid_t nsamples_did;
    hid_t visdata_fs;
    hid_t flags_fs;
    hid_t nsamples_fs;
} hdf5_id_t;  

void close_file(hdf5_id_t *id, double file_stop_t, double file_duration, uint64_t file_nblts) {
    hid_t dataset_id;
    dataset_id = H5Dopen(id->extra_keywords_gid, "stopt", H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_stop_t);
    H5Dclose(dataset_id);
    dataset_id = H5Dopen(id->extra_keywords_gid, "duration", H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_duration);
    H5Dclose(dataset_id);
    dataset_id = H5Dopen(id->header_gid, "Nblts", H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_nblts);
    H5Dclose(dataset_id);
    // Close datasets
    H5Dclose(id->visdata_did);
    H5Dclose(id->flags_did);
    H5Dclose(id->nsamples_did);
    // Close groups
    H5Gclose(id->extra_keywords_gid);
    H5Gclose(id->header_gid);
    H5Gclose(id->data_gid);
    // Close file
    H5Fflush(id->file_id, H5F_SCOPE_GLOBAL);
    H5Fclose(id->file_id);
}

/*
Create an extensible dataset for visdata, visdata_diff, flags, nsamples
see https://gist.github.com/simleb/5205083/
*/
// see https://support.hdfgroup.org/services/contributions.html
# define FILTER_H5_LZF 32000
# define FILTER_H5_BITSHUFFLE 32008
void make_extensible_hdf5(hdf5_id_t *id)
{
    hsize_t dims[N_DATA_DIMS] = {0, VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES, N_CHAN_PROCESSED, N_STOKES};
    hsize_t max_dims[N_DATA_DIMS] = {16, VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES, N_CHAN_PROCESSED, N_STOKES};
    hsize_t chunk_dims[N_DATA_DIMS] = {1, 1, N_CHAN_PROCESSED, N_STOKES};
    
    hid_t file_space = H5Screate_simple(N_DATA_DIMS, dims, max_dims);
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);

    //hid_t plist_lzf = H5Pcreate(H5P_DATASET_CREATE);
    //H5Pset_chunk(plist_lzf, N_DATA_DIMS, chunk_dims);
    //H5Pset_shuffle(plist_lzf);
    //H5Pset_filter(plist_lzf, H5PY_FILTER_LZF, H5Z_FLAG_OPTIONAL, 0, NULL);

    H5Pset_layout(plist, H5D_CHUNKED);

    H5Pset_chunk(plist, N_DATA_DIMS, chunk_dims);
    
    // Now we have the dataspace properties, create the datasets
    id->visdata_did = H5Dcreate(id->data_gid, "visdata", complex_id, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    if (id->visdata_did < 0) {
        hashpipe_error(__FUNCTION__, "Failed to create visdata dataset");
        pthread_exit(NULL);
    }

    id->nsamples_did = H5Dcreate(id->data_gid, "nsamples", H5T_STD_I64LE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    if (id->nsamples_did < 0) {
        hashpipe_error(__FUNCTION__, "Failed to create nsamples dataset");
        pthread_exit(NULL);
    }

    id->flags_did = H5Dcreate(id->data_gid, "flags", H5T_NATIVE_HBOOL, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    if (id->flags_did < 0) {
        hashpipe_error(__FUNCTION__, "Failed to create flags dataset");
        pthread_exit(NULL);
    }

    H5Pclose(plist);
    H5Sclose(file_space);
}

hid_t open_hdf5_from_template(char * sourcename, char * destname)
{
    int read_fd, write_fd;
    struct stat stat_buf;
    off_t offset = 0;
    hid_t status;
    
    read_fd = open(sourcename, O_RDONLY);
    if (read_fd < 0) {
        hashpipe_error(__FUNCTION__, "error opening %s", sourcename);
        pthread_exit(NULL);
    }

    // Get the file size
    fstat(read_fd, &stat_buf);

    write_fd = open(destname, O_WRONLY | O_CREAT, stat_buf.st_mode);
    if (write_fd < 0) {
        hashpipe_error(__FUNCTION__, "error opening %s", destname);
        pthread_exit(NULL);
    }

    sendfile(write_fd, read_fd, &offset, stat_buf.st_size);
    
    close(read_fd);
    close(write_fd);

    status = H5Fopen(destname, H5F_ACC_RDWR, H5P_DEFAULT);
    if (status < 0) {
        hashpipe_error(__FUNCTION__, "error opening %s as HDF5 file", destname);
        pthread_exit(NULL);
    }

    return status;
}


void start_file(hdf5_id_t *id, char *template_fname, char *hdf5_fname, uint64_t file_obs_id, double file_start_t) {
    hid_t dataset_id;

    id->file_id = open_hdf5_from_template(template_fname, hdf5_fname);
    // Open HDF5 header groups and create data group

    id->header_gid = H5Gopen(id->file_id, "Header", H5P_DEFAULT);
    if (id->header_gid < 0) {
        hashpipe_error(__FUNCTION__, "Failed to open Header");
        pthread_exit(NULL);
    }
    
    id->extra_keywords_gid = H5Gopen(id->header_gid, "extra_keywords", H5P_DEFAULT);
    if (id->extra_keywords_gid < 0) {
        hashpipe_error(__FUNCTION__, "Failed to open Header/extra_keywords");
        pthread_exit(NULL);
    }

    id->data_gid = H5Gcreate(id->file_id, "Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (id->data_gid < 0) {
        hashpipe_error(__FUNCTION__, "Failed to create Data group");
        pthread_exit(NULL);
    }
    // Create the extensible "Data" group datasets. This function
    // assigns all the dataset ids to the id struct
    make_extensible_hdf5(id);
    
    // Write meta-data values we know at file-open
    dataset_id = H5Dopen(id->extra_keywords_gid, "obs_id", H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_obs_id);
    H5Dclose(dataset_id);
    dataset_id = H5Dopen(id->extra_keywords_gid, "startt", H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_start_t);
    H5Dclose(dataset_id);
}

void extend_datasets(hdf5_id_t *id, int n) {
    hsize_t dims[N_DATA_DIMS] = {n, VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES, N_CHAN_PROCESSED, N_STOKES};
    H5Dset_extent(id->visdata_did, dims);
    H5Dset_extent(id->flags_did, dims);
    H5Dset_extent(id->nsamples_did, dims);
    id->visdata_fs  = H5Dget_space(id->visdata_did);
    id->flags_fs    = H5Dget_space(id->flags_did);
    id->nsamples_fs = H5Dget_space(id->nsamples_did);
}

void close_filespaces(hdf5_id_t *id) {
    H5Sclose(id->visdata_fs);
    H5Sclose(id->flags_fs);
    H5Sclose(id->nsamples_fs);
}



/*
Write an n_baselines x n_stokes data block to dataset `id` at time position `t`, channel number `c`
*/
void write_channels(hdf5_id_t *id, hsize_t t, hsize_t b, hid_t mem_space, uint64_t *visdata_buf)
{
    hsize_t start[N_DATA_DIMS] = {t, b, 0, 0};
    hsize_t count[N_DATA_DIMS] = {1, N_BL_PER_WRITE, N_CHAN_PROCESSED, N_STOKES};
    H5Sselect_hyperslab(id->visdata_fs, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dwrite(id->visdata_did, complex_id, mem_space, id->visdata_fs, H5P_DEFAULT, visdata_buf);
}

/*
Turn an mcnt into a UNIX time in double-precision.
*/
double mcnt2time(uint64_t mcnt, uint32_t sync_time)
{
    return (sync_time + mcnt) * (N_CHAN_TOTAL_GENERATED / (double)FENG_SAMPLE_RATE);
}

/*
Given an intire input buffer --
even/odd x 1 time x N_chans x N_bls x  N_stokes x 2 (real/imag)
Write --
even/odd-sum x 1 bl x N_chans x N_stokes x 2 (real/imag) into `out_sum`.
even/odd-diff dx 1 bl x N_chans x N_stokes x 2 (real/imag) into `out_diff`.
Choose the baseline with the index `b`
This function sums over CATCHER_SUM_CHANS as it transposes, and computes even/odd sum/diff.

*/
#define corr_databuf_data_idx(c, b) (2*((c*VIS_MATRIX_ENTRIES_PER_CHAN) + (N_STOKES*b)))
void transpose_bl_chan(int32_t *in, int32_t *out_sum, int32_t *out_diff, int b) {
    
    int i, chan, output_chan;
    __m256i sum_even[N_BL_PER_WRITE], sum_odd[N_BL_PER_WRITE];
    __m256i val_even[N_BL_PER_WRITE], val_odd[N_BL_PER_WRITE];
    // Wasteful, but initialize so that the compiler is happy...
    for(i=0; i<N_BL_PER_WRITE; i++) {
        sum_even[i] = _mm256_set_epi64x(0ULL,0ULL,0ULL,0ULL);
        sum_odd[i]  = _mm256_set_epi64x(0ULL,0ULL,0ULL,0ULL);
    }

    int32_t *in_even = in + corr_databuf_data_idx(0,b);
    int32_t *in_odd = in_even + 2*VIS_MATRIX_ENTRIES;
    __m256i *in_even256 = (__m256i *)in_even;
    __m256i *in_odd256  = (__m256i *)in_odd;
    __m256i *out_sum256 = (__m256i *)out_sum;
    __m256i *out_diff256 = (__m256i *)out_diff;
    int32_t temp[8];

    for (chan=0; chan<N_CHAN_TOTAL; chan++) {
        //val_odd  = _mm256_load_si256(in_odd256 + 1);
        for(i=0; i<N_BL_PER_WRITE; i++) {
            val_even[i] = _mm256_load_si256(in_even256 + i);
            val_odd[i] = _mm256_load_si256(in_odd256 + i);
            _mm256_store_si256((__m256i *)temp, val_even[i]);
            //fprintf(stdout, "%d, %d, %d, %d, %d, %d, %d, %d\n", temp[0], temp[1], temp[2], temp[3], temp[4], temp[5], temp[6], temp[7]);
        }
        if ((chan % CATCHER_CHAN_SUM) == 0) {
            // Start new accumulation (count to VIS_MATRIX_ENTRIES_PER_CHAN*2 for real/imag)
            for(i=0; i<N_BL_PER_WRITE; i++) {
                sum_even[i] = val_even[i];
                sum_odd[i] = val_odd[i];
            }
        } else {
            // Add to existing accumulation
            for(i=0; i<N_BL_PER_WRITE; i++) {
                sum_even[i] = _mm256_add_epi32(sum_even[i], val_even[i]);
                sum_odd[i] = _mm256_add_epi32(sum_odd[i], val_odd[i]);
            }
        }

        in_even256 += (VIS_MATRIX_ENTRIES_PER_CHAN >> 2);
        in_odd256  += (VIS_MATRIX_ENTRIES_PER_CHAN >> 2);
        // If this is the last channel in the sum, take the even/odd sum/diff, and write to the
        // output buffer
        if ((chan % CATCHER_CHAN_SUM) == (CATCHER_CHAN_SUM-1)) {
            output_chan = chan / CATCHER_CHAN_SUM;
            // Compute sums/diffs
            for(i=0; i<N_BL_PER_WRITE; i++) {
                _mm256_store_si256(out_sum256  + i*N_CHAN_PROCESSED + output_chan, _mm256_add_epi32(sum_even[i], sum_odd[i]));
                _mm256_store_si256(out_diff256 + i*N_CHAN_PROCESSED + output_chan, _mm256_sub_epi32(sum_even[i], sum_odd[i]));
            }
        }
    }
}


static int init(hashpipe_thread_args_t *args)
{
    //hashpipe_status_t st = args->st;
    if (register_lzf() < 0) {
        hashpipe_error(__FUNCTION__, "error registering LZF filter");
        pthread_exit(NULL);
    }
    fprintf(stdout, "Initializing Catcher disk thread\n");

    // generate the complex data type
    complex_id = H5Tcreate(H5T_COMPOUND, 8);
    H5Tinsert(complex_id, "r", 0, H5T_STD_I32LE);
    H5Tinsert(complex_id, "i", 4, H5T_STD_I32LE);
    return 0;
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
    double gps_time;
    double julian_time;

    // How many integrations to dump to a file before starting the next one
    // This is read from shared memory
    uint32_t ms_per_file;

    // Get template filename from redis
    hashpipe_status_lock_safe(&st);
    hgets(st.buf, "HDF5TPLT", 128, template_fname);
    
    // Get time that F-engines were last sync'd
    hgetu4(st.buf, "SYNCTIME", &sync_time);

    hgetu4(st.buf, "MSPERFIL", &ms_per_file);

    // Init status variables
    hputi8(st.buf, "DISKMCNT", 0);
    hashpipe_status_unlock_safe(&st);
    
    fprintf(stdout, "Catcher using header template %s\n", template_fname);
    fprintf(stdout, "Catcher using sync time %u\n", sync_time);
    fprintf(stdout, "Catcher recording %ums per file\n", ms_per_file);

    /* Loop */
    int32_t *db_in32;
    int rv;
    int curblock_in=0;
    int bl;
    double curr_file_time = -1.0;
    double file_start_t, file_stop_t, file_duration;
    int64_t file_obs_id, file_nblts=0, file_nts=0;
    float gbps, min_gbps;

    struct timespec start, finish;

    hdf5_id_t sum_file, diff_file;
    
    int32_t *bl_buf_sum  = (int32_t *)aligned_alloc(32, N_BL_PER_WRITE * N_CHAN_PROCESSED * N_STOKES * 2 * sizeof(int32_t));
    int32_t *bl_buf_diff = (int32_t *)aligned_alloc(32, N_BL_PER_WRITE * N_CHAN_PROCESSED * N_STOKES * 2 * sizeof(int32_t));

    // Allocate an array of bools for flags and n_samples
    hbool_t *flags = (hbool_t *)malloc(N_BL_PER_WRITE * N_CHAN_PROCESSED * sizeof(hbool_t));
    uint64_t *nsamples = (uint64_t *)malloc(N_BL_PER_WRITE * N_CHAN_PROCESSED* sizeof(uint64_t));

    // Define the memory space used by these buffers for HDF5 access
    // We write 1 baseline-time at a time
    hsize_t dims[N_DATA_DIMS] = {1, N_BL_PER_WRITE, N_CHAN_PROCESSED, N_STOKES};
    hid_t mem_space = H5Screate_simple(N_DATA_DIMS, dims, NULL);

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
        fprintf(stdout, "Processing new block with: mcnt: %lu\n", db_in->block[curblock_in].header.mcnt);
        fprintf(stdout, "Processing new block with: gps_time: %lf\n", gps_time);
        julian_time = 2440587.5 + (gps_time / (double)(86400.0));

        if ((curr_file_time < 0) || (1000*(gps_time - curr_file_time) > ms_per_file)) {
            fprintf(stdout, "New file trigger: gps_time: %lf, curr_file_time: %lf\n", gps_time, curr_file_time);
            // If a file is open, finish its meta-data and close it.
            if (curr_file_time >= 0) {
                fprintf(stdout, "Closing datasets and files\n");
                close_file(&sum_file, file_stop_t, file_duration, file_nblts);
                close_file(&diff_file, file_stop_t, file_duration, file_nblts);
            }

            // And now start a new file
            curr_file_time = gps_time;
            file_nblts = 0;
            file_nts = 0;
            file_start_t = gps_time;
            file_obs_id = (int64_t)gps_time;
            sprintf(hdf5_fname, "zen.%7.5lf.uvh5", julian_time);
            fprintf(stdout, "Opening new file %s\n", hdf5_fname);
            start_file(&sum_file, template_fname, hdf5_fname, file_obs_id, file_start_t);
            sprintf(hdf5_fname, "zen.%7.5lf.diff.uvh5", julian_time);
            fprintf(stdout, "Opening new file %s\n", hdf5_fname);
            start_file(&diff_file, template_fname, hdf5_fname, file_obs_id, file_start_t);
        }

        // Update time and sample counters
        file_stop_t = gps_time;
        file_duration = file_stop_t - file_start_t; //really want a +1 * acc_len here
        file_nblts += VIS_MATRIX_ENTRIES_PER_CHAN;
        file_nts += 1;

        // extend the datasets with time axes and update filespace IDs
        extend_datasets(&sum_file, file_nts);
        extend_datasets(&diff_file, file_nts);
            
        // TODO: Write lst_array, time_array, uvw_array, zenith_dec, zenith_ra

        // Sum over channels, compute even/odd sum/diffs, and get data on a per-baseline basis
        fprintf(stdout, "Writing integration\n");
        for(bl=0; bl<(VIS_MATRIX_ENTRIES_PER_CHAN/N_STOKES); bl+=N_BL_PER_WRITE) {
            transpose_bl_chan(db_in32, bl_buf_sum, bl_buf_diff, bl);
            //write data to file
            write_channels(&sum_file, file_nts-1, bl, mem_space, (uint64_t *)bl_buf_sum); 
            write_channels(&diff_file, file_nts-1, bl, mem_space, (uint64_t *)bl_buf_diff); 
            //write_baselines(nsamples_id, file_nts-1, chunk_index*CHAN_CHUNK_SIZE, fs_nsamples_id, mem_space, (uint64_t *)nsamples, H5T_STD_I64LE); 
            flags = flags;
            nsamples = nsamples;
        }

        fprintf(stdout, "Integration written\n");
        // Close the filespaces - leave the datasets open. We'll close those when the file is done
        close_filespaces(&sum_file);
        close_filespaces(&diff_file);

        clock_gettime(CLOCK_MONOTONIC, &finish);

        // Note processing time for this integration
        hashpipe_status_lock_safe(&st);

        hgetr4(st.buf, "DISKMING", &min_gbps);
        gbps = (float)(2*64L*VIS_MATRIX_ENTRIES/CATCHER_CHAN_SUM)/ELAPSED_NS(start,finish); //Gigabits / s
        hputr4(st.buf, "DISKGBPS", gbps);
        hputr4(st.buf, "DUMPMS", ELAPSED_NS(start,finish) / 1000000.0);
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

static hashpipe_thread_desc_t hera_catcher_disk_thread = {
    name: "hera_catcher_disk_thread",
    skey: "DISKSTAT",
    init: init,
    run:  run,
    ibuf_desc: {hera_catcher_input_databuf_create},
    obuf_desc: {NULL}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&hera_catcher_disk_thread);
}
