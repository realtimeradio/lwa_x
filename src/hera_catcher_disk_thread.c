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
//#include "lzf_filter.h"

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define N_DATA_DIMS (4)
#define N_CHAN_PROCESSED (N_CHAN_TOTAL / CATCHER_CHAN_SUM)
#define N_BL_PER_WRITE (32)

static hid_t complex_id;

typedef struct {
    hid_t file_id;
    hid_t header_gid;
    hid_t data_gid;
    hid_t extra_keywords_gid;
    hid_t visdata_did;
    hid_t flags_did;
    hid_t nsamples_did;
    hid_t time_array_did;
    hid_t integration_time_did;
    hid_t uvw_array_did;
    hid_t visdata_fs;
    hid_t flags_fs;
    hid_t nsamples_fs;
    hid_t time_array_fs;
    hid_t integration_time_fs;
    hid_t uvw_array_fs;
} hdf5_id_t;

static void close_file(hdf5_id_t *id, double file_stop_t, double file_duration, uint64_t file_nblts, uint64_t file_nts) {
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
    dataset_id = H5Dopen(id->header_gid, "Ntimes", H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_nts);
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
static void make_extensible_hdf5(hdf5_id_t *id)
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

/* Create the extensible header entries which have dimensions ~Nblts.
 * These are:
 * Header/uvw_array (Nblts x 3)
 * Header/time_array (Nblts)
 * Header/integration_time (Nblts)
 */
#define DIM1 1
#define DIM2 2
static void make_extensible_headers_hdf5(hdf5_id_t *id)
{
    hsize_t dims1[DIM1] = {0};
    hsize_t max_dims1[DIM1]   = {16 * VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES};
    hsize_t chunk_dims1[DIM1] = {VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES};

    hid_t file_space = H5Screate_simple(DIM1, dims1, max_dims1);
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);

    H5Pset_layout(plist, H5D_CHUNKED);

    H5Pset_chunk(plist, DIM1, chunk_dims1);

    // Now we have the dataspace properties, create the datasets
    id->time_array_did = H5Dcreate(id->header_gid, "time_array", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    if (id->time_array_did < 0) {
        hashpipe_error(__FUNCTION__, "Failed to create time_array dataset");
        pthread_exit(NULL);
    }

    id->integration_time_did = H5Dcreate(id->header_gid, "integration_time", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    if (id->integration_time_did < 0) {
        hashpipe_error(__FUNCTION__, "Failed to create integration_time dataset");
        pthread_exit(NULL);
    }

    H5Pclose(plist);
    H5Sclose(file_space);

    /* And now uvw_array, which has shape Nblts x 3 */
    hsize_t dims2[DIM2] = {0, 3};
    hsize_t max_dims2[DIM2]   = {16 * VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES, 3};
    hsize_t chunk_dims2[DIM2] = {VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES, 3};

    file_space = H5Screate_simple(DIM2, dims2, max_dims2);
    plist = H5Pcreate(H5P_DATASET_CREATE);

    H5Pset_layout(plist, H5D_CHUNKED);

    H5Pset_chunk(plist, DIM2, chunk_dims2);

    // Now we have the dataspace properties, create the datasets
    id->uvw_array_did = H5Dcreate(id->header_gid, "uvw_array", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    if (id->uvw_array_did < 0) {
        hashpipe_error(__FUNCTION__, "Failed to create uvw_array dataset");
        pthread_exit(NULL);
    }

    H5Pclose(plist);
    H5Sclose(file_space);
}

static hid_t open_hdf5_from_template(char * sourcename, char * destname)
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


static void start_file(hdf5_id_t *id, char *template_fname, char *hdf5_fname, uint64_t file_obs_id, double file_start_t) {
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
    // Create the extensible "Header/*" group datasets. This function
    // assigns all the dataset ids to the id struct
    make_extensible_headers_hdf5(id);
    
    // Write meta-data values we know at file-open
    dataset_id = H5Dopen(id->extra_keywords_gid, "obs_id", H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_obs_id);
    H5Dclose(dataset_id);
    dataset_id = H5Dopen(id->extra_keywords_gid, "startt", H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_start_t);
    H5Dclose(dataset_id);
}

static void extend_datasets(hdf5_id_t *id, int n) {
    hsize_t dims[N_DATA_DIMS] = {n, VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES, N_CHAN_PROCESSED, N_STOKES};
    H5Dset_extent(id->visdata_did, dims);
    H5Dset_extent(id->flags_did, dims);
    H5Dset_extent(id->nsamples_did, dims);
    id->visdata_fs  = H5Dget_space(id->visdata_did);
    id->flags_fs    = H5Dget_space(id->flags_did);
    id->nsamples_fs = H5Dget_space(id->nsamples_did);
}

static void extend_header_datasets(hdf5_id_t *id, int n) {
    hsize_t dims1[DIM1] = {n * VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES};
    hsize_t dims2[DIM2] = {n * VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES, 3};
    H5Dset_extent(id->time_array_did, dims1);
    H5Dset_extent(id->integration_time_did, dims1);
    H5Dset_extent(id->uvw_array_did, dims2);
    id->time_array_fs = H5Dget_space(id->time_array_did);
    id->integration_time_fs = H5Dget_space(id->integration_time_did);
    id->uvw_array_fs = H5Dget_space(id->uvw_array_did);
}

static void close_filespaces(hdf5_id_t *id) {
    H5Sclose(id->visdata_fs);
    H5Sclose(id->flags_fs);
    H5Sclose(id->nsamples_fs);
    H5Sclose(id->time_array_fs);
    H5Sclose(id->integration_time_fs);
    H5Sclose(id->uvw_array_fs);
}



/*
Write an n_baselines x n_stokes data block to dataset `id` at time position `t`, channel number `c`
*/
static void write_channels(hdf5_id_t *id, hsize_t t, hsize_t b, hid_t mem_space, uint64_t *visdata_buf)
{
    hsize_t start[N_DATA_DIMS] = {t, b, 0, 0};
    hsize_t count[N_DATA_DIMS] = {1, N_BL_PER_WRITE, N_CHAN_PROCESSED, N_STOKES};
    H5Sselect_hyperslab(id->visdata_fs, H5S_SELECT_SET, start, NULL, count, NULL);
    H5Dwrite(id->visdata_did, complex_id, mem_space, id->visdata_fs, H5P_DEFAULT, visdata_buf);
}

/*
Write Nbls entries into the integration_time and time_array arrays. Write Nbls x 3 entries into uvw_array
*/
static void write_extensible_headers(hdf5_id_t *id, hsize_t t, hid_t mem_space1, hid_t mem_space2, double *integration_time_buf, double *time_array_buf, double*uvw_array_buf)
{
    hsize_t start1[DIM1] = {t * (VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES)};
    hsize_t count1[DIM1] = {VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES};
    hsize_t start2[DIM2] = {t * (VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES), 0};
    hsize_t count2[DIM2] = {VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES, 3};
    H5Sselect_hyperslab(id->integration_time_fs, H5S_SELECT_SET, start1, NULL, count1, NULL);
    H5Dwrite(id->integration_time_did, H5T_NATIVE_DOUBLE, mem_space1, id->integration_time_fs, H5P_DEFAULT, integration_time_buf);
    H5Sselect_hyperslab(id->time_array_fs, H5S_SELECT_SET, start1, NULL, count1, NULL);
    H5Dwrite(id->time_array_did, H5T_NATIVE_DOUBLE, mem_space1, id->time_array_fs, H5P_DEFAULT, time_array_buf);
    H5Sselect_hyperslab(id->uvw_array_fs, H5S_SELECT_SET, start2, NULL, count2, NULL);
    H5Dwrite(id->uvw_array_did, H5T_NATIVE_DOUBLE, mem_space2, id->uvw_array_fs, H5P_DEFAULT, uvw_array_buf);
}

/*
Turn an mcnt into a UNIX time in double-precision.
*/
static double mcnt2time(uint64_t mcnt, uint32_t sync_time)
{
    return (sync_time + mcnt) * (N_CHAN_TOTAL_GENERATED / (double)FENG_SAMPLE_RATE);
}

static void compute_time_array(double time, double *time_buf)
{
    int i;
    for (i=0; i<(VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES); i++) {
        time_buf[i] = time;
    }
}

static void compute_integration_time_array(double integration_time, double *integration_time_buf)
{
    int i;
    for (i=0; i<(VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES); i++) {
        integration_time_buf[i] = integration_time;
    }
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
static void transpose_bl_chan(int32_t *in, int32_t *out_sum, int32_t *out_diff, int b) {
    
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

    for (chan=0; chan<N_CHAN_TOTAL; chan++) {
        //val_odd  = _mm256_load_si256(in_odd256 + 1);
        for(i=0; i<N_BL_PER_WRITE; i++) {
            val_even[i] = _mm256_load_si256(in_even256 + i);
            val_odd[i] = _mm256_load_si256(in_odd256 + i);
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
                _mm256_stream_si256(out_sum256  + i*N_CHAN_PROCESSED + output_chan, _mm256_add_epi32(sum_even[i], sum_odd[i]));
                _mm256_stream_si256(out_diff256 + i*N_CHAN_PROCESSED + output_chan, _mm256_sub_epi32(sum_even[i], sum_odd[i]));
            }
        }
    }
}


static int init(hashpipe_thread_args_t *args)
{
    //hashpipe_status_t st = args->st;
    //if (register_lzf() < 0) {
    //    hashpipe_error(__FUNCTION__, "error registering LZF filter");
    //    pthread_exit(NULL);
    //}
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

    // Timers for performance monitoring
    struct timespec t_start, w_start, w_stop; // transpose start, write start, write stop
    float gbps, min_gbps;

    struct timespec start, finish;
    uint64_t t_ns;
    uint64_t w_ns;
    uint64_t min_t_ns = 999999999;
    uint64_t min_w_ns = 999999999;
    uint64_t max_t_ns = 0;
    uint64_t max_w_ns = 0;
    uint64_t elapsed_t_ns = 0;
    uint64_t elapsed_w_ns = 0;
    float bl_t_ns = 0.0;
    float bl_w_ns = 0.0;
     
    // Buffers for file name strings
    char template_fname[128];
    char hdf5_fname[128];

    // Variables for sync time and computed gps time / JD
    uint32_t sync_time = 0;
    double gps_time;
    double julian_time;
    double integration_time;

    // How many integrations to dump to a file before starting the next one
    // This is read from shared memory
    uint32_t ms_per_file;

    // Init status variables
    hashpipe_status_lock_safe(&st);
    hputi8(st.buf, "DISKMCNT", 0);
    hashpipe_status_unlock_safe(&st);
    
    /* Loop */
    int32_t *db_in32;
    int rv;
    int curblock_in=0;
    int bl;
    double curr_file_time = -1.0;
    double file_start_t, file_stop_t, file_duration;
    int64_t file_obs_id, file_nblts=0, file_nts=0;

    hdf5_id_t sum_file, diff_file;
    
    int32_t *bl_buf_sum  = (int32_t *)aligned_alloc(32, N_BL_PER_WRITE * N_CHAN_PROCESSED * N_STOKES * 2 * sizeof(int32_t));
    int32_t *bl_buf_diff = (int32_t *)aligned_alloc(32, N_BL_PER_WRITE * N_CHAN_PROCESSED * N_STOKES * 2 * sizeof(int32_t));

    double *integration_time_buf = (double *)malloc((VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES) * sizeof(double));
    double *time_array_buf       = (double *)malloc((VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES) * sizeof(double));
    double *uvw_array_buf        = (double *)malloc((VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES) * 3 * sizeof(double));

    // Allocate an array of bools for flags and n_samples
    hbool_t *flags = (hbool_t *)malloc(N_BL_PER_WRITE * N_CHAN_PROCESSED * sizeof(hbool_t));
    uint64_t *nsamples = (uint64_t *)malloc(N_BL_PER_WRITE * N_CHAN_PROCESSED* sizeof(uint64_t));

    // Define the memory space used by these buffers for HDF5 access
    // We write 1 baseline-time at a time
    hsize_t dims[N_DATA_DIMS] = {1, N_BL_PER_WRITE, N_CHAN_PROCESSED, N_STOKES};
    hid_t mem_space = H5Screate_simple(N_DATA_DIMS, dims, NULL);
    // Memory spaces to Nblts-element header vectors
    hsize_t dims1[DIM1] = {VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES};
    hid_t mem_space1 = H5Screate_simple(DIM1, dims1, NULL);
    // Memory spaces to (Nblts x 3)-element header vectors
    hsize_t dims2[DIM2] = {VIS_MATRIX_ENTRIES_PER_CHAN / N_STOKES, 3};
    hid_t mem_space2 = H5Screate_simple(DIM2, dims2, NULL);

    while (run_threads()) {
        // Note waiting status,
        // query integrating status
        // and, if armed, start count
        hashpipe_status_lock_safe(&st);
        hputs(st.buf, status_key, "waiting");
        hashpipe_status_unlock_safe(&st);

        // Wait for new input block to be filled
        fprintf(stdout, "disk thread waiting for block %d to be filled\n", curblock_in);
        while ((rv=hera_catcher_input_databuf_busywait_filled(db_in, curblock_in)) != HASHPIPE_OK) {
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
        fprintf(stdout, "disk thread acquired block %d\n", curblock_in);

        // Get template filename from redis
        hashpipe_status_lock_safe(&st);
        hgets(st.buf, "HDF5TPLT", 128, template_fname);
        
        // Get time that F-engines were last sync'd
        hgetu4(st.buf, "SYNCTIME", &sync_time);

        hgetu4(st.buf, "MSPERFIL", &ms_per_file);

        // Get the integration time reported by the correlator
        hgetr8(st.buf, "INTSECS", &integration_time);

        // Got a new data block, update status
        hputs(st.buf, status_key, "writing");
        hputi4(st.buf, "DISKBKIN", curblock_in);
        hputu8(st.buf, "DISKMCNT", db_in->block[curblock_in].header.mcnt);
        hashpipe_status_unlock_safe(&st);
        
        db_in32 = (int32_t *)db_in->block[curblock_in].data;

        clock_gettime(CLOCK_MONOTONIC, &start);
        gps_time = mcnt2time(db_in->block[curblock_in].header.mcnt, sync_time);
        //fprintf(stdout, "Processing new block with: mcnt: %lu (gps time: %lf)\n", db_in->block[curblock_in].header.mcnt, gps_time);
        julian_time = 2440587.5 + (gps_time / (double)(86400.0));

        if ((curr_file_time < 0) || (1000*(gps_time - curr_file_time) > ms_per_file)) {
            fprintf(stdout, "New file trigger: gps_time: %lf, curr_file_time: %lf\n", gps_time, curr_file_time);
            // If a file is open, finish its meta-data and close it.
            if (curr_file_time >= 0) {
                fprintf(stdout, "Closing datasets and files\n");
                close_file(&sum_file, file_stop_t, file_duration, file_nblts, file_nts);
                close_file(&diff_file, file_stop_t, file_duration, file_nblts, file_nts);
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
        extend_header_datasets(&sum_file, file_nts);
        extend_header_datasets(&diff_file, file_nts);
            
        // Write this integration's entries for lst_array, time_array, uvw_array
        // TODO: compute uvw array
        compute_time_array(julian_time, time_array_buf);
        compute_integration_time_array(integration_time, integration_time_buf);
        write_extensible_headers(&sum_file, file_nts-1, mem_space1, mem_space2, integration_time_buf, time_array_buf, uvw_array_buf);
        write_extensible_headers(&diff_file, file_nts-1, mem_space1, mem_space2, integration_time_buf, time_array_buf, uvw_array_buf);

        // Sum over channels, compute even/odd sum/diffs, and get data on a per-baseline basis
        for(bl=0; bl<(VIS_MATRIX_ENTRIES_PER_CHAN/N_STOKES); bl+=N_BL_PER_WRITE) {
	    clock_gettime(CLOCK_MONOTONIC, &t_start);
            transpose_bl_chan(db_in32, bl_buf_sum, bl_buf_diff, bl);
            //write data to file
	    clock_gettime(CLOCK_MONOTONIC, &w_start);
            write_channels(&sum_file, file_nts-1, bl, mem_space, (uint64_t *)bl_buf_sum); 
            write_channels(&diff_file, file_nts-1, bl, mem_space, (uint64_t *)bl_buf_diff); 
	    clock_gettime(CLOCK_MONOTONIC, &w_stop);
            flags = flags;
            nsamples = nsamples;

	    t_ns = ELAPSED_NS(t_start, w_start);
	    w_ns = ELAPSED_NS(w_start, w_stop);
            elapsed_t_ns += t_ns;
            elapsed_w_ns += w_ns;
            min_t_ns = MIN(t_ns, min_t_ns);
            max_t_ns = MAX(t_ns, min_t_ns);
            min_w_ns = MIN(w_ns, min_w_ns);
            max_w_ns = MAX(w_ns, min_w_ns);
        }


        // Close the filespaces - leave the datasets open. We'll close those when the file is done
        close_filespaces(&sum_file);
        close_filespaces(&diff_file);

        clock_gettime(CLOCK_MONOTONIC, &finish);

        // Note processing time for this integration and update other stats
        bl_t_ns = (float)elapsed_t_ns / VIS_MATRIX_ENTRIES_PER_CHAN;
        bl_w_ns = (float)elapsed_w_ns / VIS_MATRIX_ENTRIES_PER_CHAN;

        hashpipe_status_lock_safe(&st);
        hputr4(st.buf, "DISKTBNS", bl_t_ns); 
        hputi8(st.buf, "DISKTMIN", min_t_ns); 
        hputi8(st.buf, "DISKTMAX", max_t_ns); 
        hputr4(st.buf, "DISKWBNS", bl_w_ns); 
        hputi8(st.buf, "DISKWMIN", min_w_ns); 
        hputi8(st.buf, "DISKWMAX", max_w_ns); 

        hgetr4(st.buf, "DISKMING", &min_gbps);
        gbps = (float)(2*64L*VIS_MATRIX_ENTRIES/CATCHER_CHAN_SUM)/ELAPSED_NS(start,finish); //Gigabits / s
        hputr4(st.buf, "DISKGBPS", gbps);
        hputr4(st.buf, "DUMPMS", ELAPSED_NS(start,finish) / 1000000.0);
        if(min_gbps == 0 || gbps < min_gbps) {
          hputr4(st.buf, "DISKMING", gbps);
        }
        hashpipe_status_unlock_safe(&st);

        fprintf(stdout, "disk thread freeing buf %d \n", curblock_in);
        // Mark input block as free and advance
        if(hera_catcher_input_databuf_set_free(db_in, curblock_in) != HASHPIPE_OK) {
            hashpipe_error(__FUNCTION__, "error marking databuf %d free", curblock_in);
            pthread_exit(NULL);
        }
        curblock_in = (curblock_in + 1) % CATCHER_N_BLOCKS;

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
