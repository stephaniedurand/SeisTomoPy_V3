#!/bin/csh

#------------------------------------------------
  set filename1 = $1
  set lmax = $2
  set dep1 = $3
  set model1 = $4
  set para1 = $5	
#------------------------------------------------

cd bin
echo ""
echo "COMPUTE CROSS-CORRELATION" 
echo ""
./make_corr_fromfile_bis.exe  << !
$filename1 
$lmax 
$model1 $para1 $dep1
!
echo "CROSS-CORRELATION DONE" 
echo ""

rm ../input_files/*

cd ..

