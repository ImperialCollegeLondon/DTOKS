wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/zlib-1.2.8.tar.gz
tar -xf zlib-1.2.8.tar.gz && cd zlib-1.2.8
ZDIR=/usr/local/zlib-1.2.8
./configure --prefix="/usr/local/zlib-1.2.8"
#sudo make check install
sudo make install
cd ..

wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/hdf5-1.8.13.tar.gz
tar -xf hdf5-1.8.13.tar.gz && cd hdf5-1.8.13
HDF5_DIR=/usr/local/hdf5-1.8.13
./configure --enable-shared --with-zlib=${ZDIR} --enable-hl --prefix="/usr/local/hdf5-1.8.13"
make # 2 for number of procs to be used
#sudo make check
sudo make install
cd ..

wget http://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-4.1.3.tar.gz
tar -xf netcdf-4.1.3.tar.gz && cd netcdf-4.1.3
CPPFLAGS=-I$HDF5_DIR/include LDFLAGS=-L$HDF5_DIR/lib ./configure --enable-netcdf-4 --enable-shared --enable-dap --prefix="/usr/local/netcdf-4.1.3"
# make check
make 
sudo make install
cd ..


