cd examples

rest="-L/mnt/UbuntuData2/pythiaDevelopment/lib -lpythia8 -I/mnt/UbuntuData2/pythiaDevelopment/include -std=c++20 -Wno-misleading-indentation -Wno-deprecated-declarations"

echo -e ">> running main01"
g++ main01.cc -o ../bin/main01 -O2 ${rest}
time ../bin/main01 > ../results/main01.txt 2> ../results/main01perf.txt

echo -e "\n>> running main02"
g++ main02.cc -o ../bin/main02 -O2 ${rest}
time ../bin/main02 > ../results/main02.txt 2> ../results/main02perf.txt

echo -e "\n>> running main03"
g++ main03.cc -o ../bin/main03 -O2 ${rest}
time ../bin/main03 > ../results/main03.txt 2> ../results/main03perf.txt

echo -e "\n>> running main04"
g++ main04.cc -o ../bin/main04 -O2 ${rest}
time ../bin/main04 > ../results/main04.txt 2> ../results/main04perf.txt

echo -e "\n>> running main24"
g++ main24.cc -o ../bin/main24 -O2 ${rest}
time ../bin/main24 > ../results/main24.txt 2> ../results/main24perf.txt

echo -e "\n>> running mainX01"
g++ mainX01.cc -o ../bin/mainX01 -O2 ${rest}
time ../bin/mainX01 > ../results/mainX01.txt 2> ../results/mainX01perf.txt

echo -e "\nall done"
