cd examples

rest="-L/mnt/UbuntuData2/pythiaDevelopment/lib -lpythia8 -I/mnt/UbuntuData2/pythiaDevelopment/include -std=c++20 -Wno-misleading-indentation -Wno-deprecated-declarations"

g++ main01.cc -o ../bin/main01 -O2 ${rest}
../bin/main01 > ../results/main01.txt 2> ../results/main01perf.txt

g++ main02.cc -o ../bin/main02 -O2 ${rest}
../bin/main02 > ../results/main02.txt 2> ../results/main02perf.txt

g++ main03.cc -o ../bin/main03 -O2 ${rest}
../bin/main03 > ../results/main03.txt 2> ../results/main03perf.txt

g++ main04.cc -o ../bin/main04 -O2 ${rest}
../bin/main04 > ../results/main04.txt 2> ../results/main04perf.txt

g++ main24.cc -o ../bin/main24 -O2 ${rest}
../bin/main24 > ../results/main24.txt 2> ../results/main24perf.txt
