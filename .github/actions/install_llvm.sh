#!/bin/bash
wget https://apt.llvm.org/llvm.sh
chmod u+x ./llvm.sh
yes | sudo  ./llvm.sh $1 all
rm -f ./llvm.sh