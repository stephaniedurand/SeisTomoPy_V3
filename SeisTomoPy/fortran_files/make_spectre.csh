#!/bin/csh

#------------------------------------------------
  set lmax = $1	
  set model = $2
  set para = $3
  set dep = $4
#------------------------------------------------

cd bin
echo ""
echo "COMPUTE SPECTRE" 
./make_spectre.exe  << !
$lmax $model $para $dep 
!
echo "SPECTRE DONE" 
echo ""

rm ../input_files/map_NEW*
rm ../output_files_map/map_NEW*
rm ../output_files_map/output_map*

cd ..

