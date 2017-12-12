#!/bin/csh

#------------------------------------------------
  set lmax = $1	
  set dep1 = $2
  set model1 = $3
  set para1 = $4
  set dep2 = $5
  set model2 = $6
  set para2 = $7
#------------------------------------------------

cd bin
echo ""
echo "COMPUTE CROSS-CORRELATION" 
echo ""
./make_corr.exe  << !
$lmax 
$model1 $para1 $dep1
$model2 $para2 $dep2
!
echo "CROSS-CORRELATION DONE" 
echo ""

rm ../input_files/map_NEW*
rm ../output_files_map/map_NEW*
rm ../output_files_map/output_map*

cd ..

