#!/bin/csh

rm bin/fort.111 

#------------------------------------------------
  set depth = $1	
  set model = $2
  set lmax = $3
#------------------------------------------------

cd bin
echo ""
echo "COMPUTE MAP" 
./make_map.exe << !
$depth $model $lmax
!
echo "MAP DONE" 
echo ""

cd ..


