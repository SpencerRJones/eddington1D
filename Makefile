EXEC1 = run_eddington.exe

FC = ifx
FFLAGS =  -g -O2
CRTM_INC = crtm_v2.2.3/include
CRTM_LIB = crtm_v2.2.3/lib
INCFLAGS = -I$(CRTM_INC)
LIBFLAGS = -lcrtm -L$(CRTM_LIB)
FLAGS = $(INCFLAGS) $(LIBFLAGS)

####################################################

SOURCECODE = run_eddington.f90

OBJS =  defs_eddington.o \
	linpak.o \
	watoptic.o \
	iceoptic.o \
	absorb-clr.o \
	monortm-lut.o \
	mie.o       \
	dsd.o       \
	Fastem.o    \
	radtran.o   \
	eddington.o \
	run_eddington.o 

####################################################

$(EXEC1):	$(OBJS)
		$(FC) $(FFLAGS) $(OBJS) $(FLAGS) -o $(EXEC1)
		rm -f *.o *.mod


tags:   $(SOURCECODE)
	ctags $(SOURCECODE)

clean:	
	rm -f *.o *.mod

.SUFFIXES: .f95 .f90 .f .o
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FLAGS) $<


