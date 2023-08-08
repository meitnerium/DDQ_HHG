FC=ifort
MKL="/opt/intel/mkl/10.0.3.020/"
FOPTS=-g

wp_DDQ_HHG_v01: wp_DDQ_HHG_v01.f90 morse.o potspec_maj.o spo_maj.o cmplxmath.o
	$(FC) $(FOPTS) cmplxmath.o morse.o potspec_maj.o spo_maj.o  wp_DDQ_HHG_v01.f90 -o wp_DDQ_HHG_v01

travail: travail.f90 morse.o potspec_maj.o spo_maj.o cmplxmath.o 
	$(FC) $(FOPTS) cmplxmath.o morse.o potspec_maj.o spo_maj.o travail.f90 -o travail	

morse.o: morse.f
	$(FC) $(FOPTS) -c morse.f

potspec_maj.o:potspec_maj.f
	$(FC) $(FOPTS) -c potspec_maj.f

spo_maj.o :spo_maj.f90 cmplxmath.o 
	$(FC) $(FOPTS) cmplxmath.o -c spo_maj.f90

cmplxmath.o:cmplxmath.f
	$(FC) $(FOPTS) -c cmplxmath.f

run: travail
	./travail

fort.10: travail
	./travail

fort.11: travail
	./travail
clean:
	rm travail *.o


