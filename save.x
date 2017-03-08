#! /bin/csh -ef
#
#  Backup da tese de mestrado na luxor
#
unalias rm
rm -f scf.x intvalues.dat 
cd ..
rm -f scf-program.tar.gz scf-program.tar
tar -cvf scf-program.tar ./scf
gzip scf-program.tar
scp scf-program.tar.gz lmartinez@143.106.51.146:./
