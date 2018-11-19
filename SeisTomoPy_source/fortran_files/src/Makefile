########################################################################################################
# Intel fortran Compilateur
# OPTION A EVITER AVEC IFORT: -fast
F90      = gfortran
# Options de compilation [mode normal] 
F90FLAGS  =  
#F90FLAGS  = -check all 
#F90FLAGS  = -check all -zero -autodouble -real-size 64 -mcmodel=large -shared-intel
#F90FLAGS  = -check all -r8 -mcmodel=large -shared-intel
#F90FLAGS = -zero -autodouble -i8 -r8 -O2 -fPIC -ipo -inline all -implicitnone -fp-model precise
#F90FLAGS = -zero -autodouble -i8 -r8 -O2 -fPIC -ipo -inline all -openmp -implicitnone -fp-model precise
#F90FLAGS = -check all -zero -autodouble -i8 -r8 -O2 -fPIC -ipo -inline all -implicitnone -fp-model precise
# Options de compilation [mode paranoiaque]
#F90FLAGS = -r8 -i8 -autodouble -O0 -fno-alias -fno-fnalias -check all -implicitnone -warn all -fp-stack-check  \
           -fPIC -zero -ftrapuv -heap-arrays -gen-interface -warn truncated_source -warn argument_checking \
           -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage  
#LDFLAGS  = -fPIC
#LDFLAGS  = -openmp -fPIC
########################################################################################################
# GNU fortran compileur
#F90      = gfortran
# Options de compilation [mode normal] 
#F90FLAGS =   
#F90FLAGS =  -mcmodel=large 
# Options de compilation [mode paranoiaque]
#F90FLAGS = -fbounds-check
#-finit-local-zero -fbounds-check -fimplicit-none \
           -pedantic -pedantic-errors -Wall -Waliasing -Wampersand -Warray-bounds  \
           -Wcharacter-truncation -Wconversion -Wline-truncation -Wintrinsics-std  \
           -Wsurprising -Wno-tabs -Wunderflow -Wintrinsic-shadow \
           -Wno-align-commons -Werror -Wunused-parameter -fdefault-real-8
#LDFLAGS  =  -fPIC
########################################################################################################
# Libraries to Link
LIB     =
LIBDIR  = 
#
MKLROOT    = 
MKLPATH    =
MKLINCLUDE =
#
MKL     = 
########################################################################################################
O = ../obj
S = ../src
B = ../bin
#
SRC      = $S/spectral.f $S/make_spectre_fromfile.f $S/make_spectre.f $S/make_corr_fromfile.f $S/make_corr_fromfile_bis.f $S/make_corr.f $S/subroutine_HS_3D2016.f $S/subroutine_HS_YOURMODEL.f $S/subroutine_HS_SP12RTS.f $S/subroutine_HS_SEMUCB.f $S/subroutine_HS_S362WMANIM.f $S/subroutine_HS_S40RTS.f $S/subroutine_HS_SGLOBE.f $S/subroutine_HS_cross_180.f $S/cross_180.f $S/subroutine_spline.f $S/subroutine_HS_map_3D2016.f $S/subroutine_HS_map_YOURMODEL.f $S/subroutine_HS_map_SP12RTS.f $S/subroutine_HS_map_SEMUCB.f $S/subroutine_HS_map_S362WMANIM.f $S/subroutine_HS_map_S40RTS.f $S/subroutine_HS_map_SGLOBE.f $S/subroutine_HS_map.f $S/make_map.f
OBJ      = $O/spectral.o $O/make_spectre_fromfile.o $O/make_spectre.o $O/make_corr_fromfile.o $O/make_corr.o $O/make_corr_fromfile_bis.o $O/subroutine_HS_3D2016.o $O/subroutine_HS_YOURMODEL.o $O/subroutine_HS_SP12RTS.o $O/subroutine_HS_SEMUCB.o $O/subroutine_HS_S362WMANIM.o $O/subroutine_HS_S40RTS.o $O/subroutine_HS_SGLOBE.o $O/subroutine_HS_cross_180.o $O/cross_180.o $O/subroutine_spline.o $O/subroutine_HS_map_3D2016.o $O/subroutine_HS_map_YOURMODEL.o $O/subroutine_HS_map_SP12RTS.o $O/subroutine_HS_map_SEMUCB.o $O/subroutine_HS_map_S362WMANIM.o $O/subroutine_HS_map_S40RTS.o $O/subroutine_HS_map_SGLOBE.o $O/subroutine_HS_map.o $O/make_map.o


LIB_corr =  $O/spectral.o \
            $O/make_corr.o 

LIB_corr_fromfile_bis =  $O/spectral.o \
               $O/make_corr_fromfile_bis.o

LIB_corr_fromfile =  $O/spectral.o \
               $O/make_corr_fromfile.o

LIB_spectre_fromfile =  $O/spectral.o \
               $O/make_spectre_fromfile.o

LIB_spectre =  $O/spectral.o \
               $O/make_spectre.o 

LIB_cross_180 =  $O/subroutine_HS_S40RTS.o \
                 $O/subroutine_HS_S362WMANIM.o \
                 $O/subroutine_HS_SEMUCB.o \
                 $O/subroutine_HS_YOURMODEL.o \
                 $O/subroutine_HS_SGLOBE.o \
                 $O/subroutine_HS_SP12RTS.o \
                 $O/subroutine_HS_3D2016.o \
                 $O/subroutine_HS_cross_180.o \
                 $O/cross_180.o 

LIB_map =  $O/subroutine_HS_map_SEMUCB.o \
           $O/subroutine_HS_map_S362WMANIM.o \
           $O/subroutine_HS_map_S40RTS.o \
           $O/subroutine_HS_map_YOURMODEL.o \
           $O/subroutine_HS_map_SP12RTS.o \
           $O/subroutine_HS_map_3D2016.o \
           $O/subroutine_HS_map_SGLOBE.o \
           $O/subroutine_spline.o \
           $O/subroutine_HS_map.o \
           $O/make_map.o 

$O/%.o: $S/%.f
	$(F90) -o $@ -c $< $(F90FLAGS)

$B/cross_180.exe: $(LIB_cross_180)
	$(F90) -o $@ $^ $(LDFLAGS) $(LIBDIR) $(LIB) $(MKL)

$B/make_map.exe: $(LIB_map)
	$(F90) -o $@ $^ $(LDFLAGS) $(LIBDIR) $(LIB) $(MKL)

$B/make_spectre.exe: $(LIB_spectre)
	$(F90) -o $@ $^ $(LDFLAGS) $(LIBDIR) $(LIB) $(MKL)

$B/make_spectre_fromfile.exe: $(LIB_spectre_fromfile)
	$(F90) -o $@ $^ $(LDFLAGS) $(LIBDIR) $(LIB) $(MKL)

$B/make_corr.exe: $(LIB_corr)
	$(F90) -o $@ $^ $(LDFLAGS) $(LIBDIR) $(LIB) $(MKL)

$B/make_corr_fromfile.exe: $(LIB_corr_fromfile)
	$(F90) -o $@ $^ $(LDFLAGS) $(LIBDIR) $(LIB) $(MKL)

$B/make_corr_fromfile_bis.exe: $(LIB_corr_fromfile_bis)
	$(F90) -o $@ $^ $(LDFLAGS) $(LIBDIR) $(LIB) $(MKL)

########################################################################################################

all: 	$B/cross_180.exe \
	$B/make_map.exe \
	$B/make_spectre.exe \
	$B/make_spectre_fromfile.exe \
	$B/make_corr.exe \
	$B/make_corr_fromfile.exe \
	$B/make_corr_fromfile_bis.exe

cleanall: clean
	rm -f ../bin/* ../obj/* ../output_files_corr/* ../output_files_map/* ../output_files_cross/* ../output_files_spectre/* ../input_files/* ../*.xy ../*.xyz ../*.res

clean:
	rm -f ../obj/* ../bin/* ../output_files_corr/* ../output_files_map/* ../output_files_cross/* ../output_files_spectre/* ../input_files/* ../*.xy ../*.xyz ../*.res

