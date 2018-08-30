#!/bin/bash
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
cd build
make
cd ..
