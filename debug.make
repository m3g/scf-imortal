#
#
#  ############################################################
#  ##                                                        ##
#  ##  compile.make  --  compile each of the TINKER modules  ##
#  ##                 (Linux/GNU g77 Version)                ##
#  ##                                                        ##
#  ############################################################
#
mv -f fixpoint.f ./memopt/fixpoint-last.f
mv -f integrals.f ./memopt/integrals-last.f
if( $1 == "memory" ) then
  cp ./memopt/fixpoint-memory.f ./fixpoint.f
  cp ./memopt/integrals-memory.f ./integrals.f
else if ( $1 == "multiple" ) then
  cp ./memopt/fixpoint-multiplefiles.f ./fixpoint.f
  cp ./memopt/integrals-multiplefiles.f ./integrals.f
else
  cp ./memopt/fixpoint-onefile.f ./fixpoint.f
  cp ./memopt/integrals-onefile.f ./integrals.f
endif
#
# Compiling:
#
g77 -c -O3 -ffast-math main.f
g77 -c -O3 -ffast-math initial.f
g77 -c -O3 -ffast-math getinp.f
g77 -c -O3 -ffast-math readxyz.f
g77 -c -O3 -ffast-math readbase.f
g77 -c -O3 -ffast-math integrals.f
g77 -c -O3 -ffast-math fixpoint.f
g77 -c -O3 -ffast-math fpacc.f
g77 -c -O3 -ffast-math auxiliar.f
g77 -c -O3 -ffast-math fmm-accelerated.f 
g77 -c -O3 -ffast-math fmm-debug.f 
#g77 -c -O3 -ffast-math fmm-simple.f 
#
# Creating object library:
# 
ar -crusv libtinker.a \
main.o \
initial.o \
getinp.o \
readxyz.o \
readbase.o \
integrals.o \
fixpoint.o \
fpacc.o \
auxiliar.o \
fmm-accelerated.o \
fmm-debug.o 
#fmm-simple.o 
#
# Linking
#
g77  -o scf.x main.o libtinker.a
rm -f *.o libtinker.a
