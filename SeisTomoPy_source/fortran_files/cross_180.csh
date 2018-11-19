#!/bin/csh

#------------------------------------------------
  set clat = $1	
  set clon = $2
  set az = $3
  set dep = $4
  set width = $5
  set model = $6
  set lmax  = $7
#------------------------------------------------

cd bin
echo ""
echo "COMPUTE CROSS-SECTION" 
./cross_180.exe  << !
$clat $clon $az $dep $width $model $lmax
!
echo "CROSS-SECTION DONE" 
echo ""

cd ..

