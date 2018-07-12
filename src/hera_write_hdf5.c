#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <hdf5.h>

#include "hera_hdf5.h"

#define FILENAME "/tmp/foo.h5"

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

void foo_write_hdf5()
{
  hid_t file_id, dataspace_id, dataset_id;
  hsize_t dims[4];
  herr_t status;
  
  file_id = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dims[0] = 1;
  dims[1] = 2;
  dims[2] = 3;
  dims[3] = 4;
  dataspace_id = H5Screate_simple(4, dims, NULL);
  
  dataset_id = H5Dcreate2(file_id, "bar", H5T_STD_I32BE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);
  status = H5Fclose(file_id);
  fprintf(stderr, "%d", (int)status);
  
}

int main(int argc, char *argv[])
{
  foo_write_hdf5();
  return 0;
}
