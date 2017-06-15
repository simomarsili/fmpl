PROG = mpl
FC = gfortran
FFLAGS = -O3 -fno-automatic
#FFLAGS = -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
#FFLAGS = -O3 
LDFLAGS = 
SRCS = kinds.f90 dvmlm.f constants.f90 parser.f90 command_line.f90 data.f90 model.f90 dvmlm_wrapper.f90 scrs.f90 mpl.f90 
OBJS = kinds.o dvmlm.o constants.o parser.o command_line.o data.o model.o dvmlm_wrapper.o scrs.o mpl.o 

$(PROG): $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.SUFFIXES: $(SUFFIXES) .f90 .f

.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

%.o: %.mod

kinds.o: kinds.f90
dvmlm.o: dvmlm.f
constants.o: constants.f90 kinds.o
parser.o: parser.f90
scrs.o: scrs.f90 kinds.o
command_line.o: command_line.f90 kinds.o constants.o
data.o: data.f90 kinds.o constants.o kinds.o parser.o
model.o: model.f90 kinds.o data.o 
dvmlm_wrapper.o: dvmlm_wrapper.f90 dvmlm.o model.o data.o kinds.o
mpl.o: mpl.f90 kinds.o constants.o parser.o command_line.o data.o model.o scrs.o dvmlm_wrapper.o 

.PHONY: clean realclean debug showtargets

# clean: remove useless files, but keep executables
clean:
	$(RM) core TAGS ?*[~#] *.o __* ...* *.mod

# realclean: remove all regenerable files, including executables
realclean: clean
	$(RM) $(PROG) $(OBJS)

# debug: debug options (enable D comments lines, disable -O for ABSOFT)
debug:
	$(MAKE) "FFLAGS=-g -fbounds-check" "OSTYPE=$(OSTYPE:-gnu=)"

# showtargets: list the most important targets of this makefile
showtargets:
	@ echo clean realclean debug fcheck fdepend 

