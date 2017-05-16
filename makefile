FC=gfortran
#FC=ifort
#FFLAG1= -mcmodel=medium -ffixed-line-length-none -fbounds-check
FFLAG1=-lfftw3 -I /usr/include -L /usr/lib64
objects1=preparation.o clogc.o filter.o taper.o do_norm.o smooth.o smoothf.o taperf.o do_whiten.o sacio.o
objects2=pre_processing.o clogc.o filter.o taper.o do_norm.o smooth.o smoothf.o taperf.o do_whiten.o sacio.o
objects3=pre_processing_2017_01_06.o clogc.o filter.o taper.o do_norm.o smooth.o smoothf.o taperf.o do_whiten_2017_01_06.o sacio.o
objects4=pre_processing_2017_03_08.o clogc.o filter_2017_03_08.o taper.o do_norm_2017_03_08.o smooth_2017_03_08.o smoothf.o taperf.o do_whiten_2017_03_08.o sacio.o
objects5=pre_processing_2017_04_11.o clogc.o filter_2017_03_08.o taper.o do_norm_2017_03_08.o smooth_2017_03_08.o smoothf.o taperf.o do_whiten_2017_03_08.o sacio.o
objects6=pre_processing_2017_04_19.o clogc.o filter_2017_04_19.o taper.o do_norm_2017_04_19.o smooth_2017_03_08.o smoothf.o taperf.o do_whiten_2017_04_19.o sacio.o
executable=preparation
all:sacio.mod preparation pre_processing  \
	pre_processing_2017_01_06 pre_processing_2017_03_08 \
	pre_processing_2017_04_11 pre_processing_2017_04_19 
.f.o:
	$(FC) $(FFLAG1) $< -c
%.o:%.f90
	$(FC) $(FFLAG1) $< -c
sacio.mod:sacio.f90
	$(FC) $< -c
preparation:$(objects1)
	$(FC) $(FFLAG1) $(objects1) -o $@
pre_processing:$(objects2)
	$(FC) $(FFLAG1) $(objects2) -o $@
pre_processing_2017_01_06:$(objects3)
	$(FC) $(FFLAG1) $(objects3) -o $@
pre_processing_2017_03_08:$(objects4)
	$(FC) $(FFLAG1) $^ -o $@
pre_processing_2017_04_11:$(objects5)
	$(FC) $(FFLAG1) $^ -o $@
pre_processing_2017_04_19:$(objects6)
	$(FC) $(FFLAG1) $^ -o $@
install:
	-cp preparation pre_processing pre_processing_2017_01_06 ../../bin
	cp pre_processing_2017_03_08 pre_processing_2017_04_11 ../../bin/
	cp pre_processing_2017_04_19 ../../bin
uninstall:
	-rm ../bin/pre_processing ../bin/preparation
clean:
	-rm *.o *.mod 
