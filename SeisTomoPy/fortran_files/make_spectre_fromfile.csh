#!/bin/csh

#------------------------------------------------
  set filename = $1
  set lmax = $2	
  echo ${filename}
#------------------------------------------------

cd bin
echo ""
echo "COMPUTE SPECTRE FROM FILE" 
echo ""
./make_spectre_fromfile.exe  << !
${filename}
$lmax 
!
echo "SPECTRE FROM FILE DONE" 
echo ""

#rm ../input_files/*

cd ..

