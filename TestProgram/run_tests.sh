#!/bin/bash

../Make/Utils/write_mkflags.pl -f ../Make/MkFlags $@

echo Cleaning...
( cd .. && make cleanall )

echo Building...
make -j2 tests

echo Run Tests...
make runalltests