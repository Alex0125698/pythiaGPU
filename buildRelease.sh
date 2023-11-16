clear
./configure --enable-debug --cxx-common="-O2 -std=c++20 -Wno-misleading-indentation -Wno-deprecated-declarations"
make clean
make -j16


# WARNING: pythia8-config doesn't actually work
# g++ main01.cc -o main01 `../bin/pythia8-config --cppflags --libs`

# just run example using
cd examples
g++ main01.cc -o main01 -O2 -L/mnt/UbuntuData2/pythiaGPU/lib -lpythia8 -I/mnt/UbuntuData2/pythiaGPU/include -std=c++20 -Wno-misleading-indentation -Wno-deprecated-declarations
./main01 > main01.txt
