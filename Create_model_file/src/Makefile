########################################################################################################
# Intel fortran Compilateur
# OPTION A EVITER AVEC IFORT: -fast
#F90      = ifort
# Options de compilation [mode normal] 
F90FLAGS  = -check all
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
F90      = gfortran
# Options de compilation [mode normal] 
F90FLAGS =  
#F90FLAGS =  -fbounds-check
# Options de compilation [mode paranoiaque]
#F90FLAGS = -finit-local-zero -fbounds-check -fimplicit-none \
           -pedantic -pedantic-errors -Wall -Waliasing -Wampersand -Warray-bounds  \
           -Wcharacter-truncation -Wconversion -Wline-truncation -Wintrinsics-std  \
           -Wsurprising -Wno-tabs -Wunderflow -Wintrinsic-shadow \
           -Wno-align-commons -Werror -Wunused-parameter -fdefault-real-8
#LDFLAGS  =  -fPIC
########################################################################################################
# Libraries to Link
LIB     = 
LIBDIR  =
########################################################################################################
#
SRC      = spectral.f create_model_file.f
OBJ      = $(SRC:.f=.o)

main: $(OBJ)
	$(F90) -o $@ $^ $(LDFLAGS) $(LIBDIR) $(LIB)

%.o: %.f
	$(F90) -o $@ -c $< $(F90FLAGS)

cleanall: clean
	rm -f main *.o *.mod

clean:
	rm -f *.o *.mod main
