#!/bin/csh

#------------------------------------------------
  set filename1 = $1
  set filename2 = $2  
  set lmax = $3	
#------------------------------------------------

cd bin
echo ""
echo "COMPUTE CROSS-CORRELATION" 
echo ""
./make_corr_fromfile.exe  << !
$filename1 
$filename2
$lmax 
!
echo "CROSS-CORRELATION DONE" 
echo ""

cd ..

