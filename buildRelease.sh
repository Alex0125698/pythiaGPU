clear
./configure --enable-debug --cxx-common="-O2 -std=c++20 -Wno-misleading-indentation -Wno-deprecated-declarations"
make clean
make -j16
