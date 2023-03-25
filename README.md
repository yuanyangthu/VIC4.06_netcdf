# VIC4.06_netcdf
Variable Infiltration Capacity Model (VIC) version 4.06 with NETCDF I/O. 

The source code is in /src, and the example global parameter file is "global_uw_0.0625d.txt". The users need to change the locations of the other input/output files and sets parameters that govern the simulation (e.g., start/end dates, modes of operation)

Compiling

cd ./src

make clean

make

Run VIC

./src/VIC_dev.exe -g global_uw_0.0625.txt (global parameter file)


