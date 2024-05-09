# This makefile creates object files for ODRPACK95, and
# also runs each of the test problems.


.SUFFIXES: .f .o .a .out .mod
F95 ?= f95 # specify compiler name as appropriate
F95OPT = # specify desired compiler options here
L = # specify options for library files and/or directories

# Specify what objects are to be created.
# OBJECTS = The object files to be created
#         = odr.o lpkbls.o
OBJECTS = odr.o lpkbls.o



# Test installation...

tests: TESTS
TESTS: drive1.out drive2.out drive3.out drive4.out test.out


# Create ODRPACK95 objects...

real_precision.o: real_precision.f
	$(F95) $(F95OPT) -c $<

lpkbls.o: lpkbls.f real_precision.o
	$(F95) $(F95OPT) -c lpkbls.f

odr.o: odr.f real_precision.o
	$(F95) $(F95OPT) -c odr.f

guide.ps: guide.dvi
	dvips -t letter -o guide.ps guide.dvi

guide.dvi: guide.tex
	tex guide.tex


# Run test problems...

drive1.out: drive1.f $(OBJECTS) data1.dat
	cp data1.dat DATA1
	$(F95) $(F95OPT) drive1.f $(OBJECTS) $L; ./a.out
	mv REPORT1 $@; rm -f DATA1 drive1.o a.out

drive2.out: drive2.f $(OBJECTS) data2.dat
	cp data2.dat DATA2
	$(F95) $(F95OPT) drive2.f $(OBJECTS) $L; ./a.out
	mv REPORT2 $@; rm -f DATA2 drive2.o a.out

drive3.out: drive3.f $(OBJECTS) data3.dat
	cp data3.dat DATA3
	$(F95) $(F95OPT) drive3.f $(OBJECTS) $L; ./a.out
	mv REPORT3 $@; rm -f DATA3 drive3.o a.out

drive4.out: drive4.f $(OBJECTS)
	$(F95) $(F95OPT) drive4.f $(OBJECTS) $L; ./a.out
	mv REPORT4 $@; rm -f drive4.o a.out

test.out: test.f $(OBJECTS)
	$(F95) $(F95OPT) test.f $(OBJECTS) $L; ./a.out
	mv REPORT $@; cat SUMMARY >> $@; rm -f test.o a.out SUMMARY

tester.out: tester.f $(OBJECTS)
	$(F95) $(F95OPT) tester.f $(OBJECTS) $L; ./a.out
	mv REPORT $@; cat SUMMARY >> $@; rm -f tester.o a.out SUMMARY

simple_example: simple_example.f90 $(OBJECTS)
	$(F95) $(F95OPT) -o simple_example simple_example.f90 $(OBJECTS) $L

restart_example: restart_example.f $(OBJECTS)
	$(F95) $(F95OPT) -o restart_example restart_example.f $(OBJECTS) $L

clean:
	rm -f *.o *.mod *.out *.dvi *.log
