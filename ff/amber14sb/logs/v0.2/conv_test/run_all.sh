#!/bin/bash

# Checks whether the conversion with qtools is ok, 
# relative to the working implementation (also checks the bonding energies)
# 
# Checks whether the working implementation is still working by comparing 
# Amber and Q single point energies of a system comprised of chain of aminoacids 
# (that are in the library) in vacuo.
#


echo "######## Qtools ########"
cd qtools
./run.sh
cd ..

echo
echo "######## AMBER ########"
cd Amber
./run.sh
cd ..

echo
echo "######## Q ########"
cd Q
./run.sh
cd ..


