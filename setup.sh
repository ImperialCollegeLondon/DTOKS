# Install dependenceis: ZLIB, HDF5, NETCDF
mkdir Dependencies
mkdir Data
cd Dependencies

. install_netcdf.sh

# Get config4cpp for configuration
git clone https://github.com/config4star/config4cpp.git

#build config 4 cpp
cd config4cpp
make
cd ..

