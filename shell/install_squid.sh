#!/bin/bash

rm -r lib
rm -r include
mkdir lib
mkdir include

mkdir include/rom_BB
mkdir include/rom_std
mkdir include/rom_sys
mkdir include/rom_ecm

cd rom_BB
make -f Makefile_squid clean
make -f Makefile_squid
cp *.h ../include/rom_BB
cp *.a ../lib/
cd ..

cd rom_std
make -f Makefile_squid clean
make -f Makefile_squid
cp *.h ../include/rom_std
cp *.a ../lib/
cd ..

cd rom_sys
make -f Makefile_squid clean
make -f Makefile_squid
cp *.h ../include/rom_sys
cp *.a ../lib/
cd ..

cd rom_ecm
make -f Makefile_squid clean
make -f Makefile_squid
cp *.h ../include/rom_ecm
cp *.a ../lib/
cd ..
