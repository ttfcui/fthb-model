
FC = gfortran
FFLAGS = -g -w -O3 -fmax-stack-var-size=65535 -ffree-line-length-none -freal-4-real-8 -fmax-errors=1

SRCS = $(patsubst %.f90, %.o, $(wildcard lifeCycle*.f90) share.f90 nrutil.f90 nrtype.f90)


lifecycle.out: $(SRCS)
	$(FC) $(FFLAGS) -fopenmp -o $@ $^

lifeCycle_TC09292016.o: lifeCycle_TC09292016.f90 lifeCycle_transition.o
	$(FC) -c $(FFLAGS) -fopenmp -o $@ $<

lifeCycle_simulate.o lifeCycle_simulate.mod: lifeCycle_simulate.f90 lifeCycle_solveDP.o lifeCycle_algs.o share.f90
	$(FC) -c $(FFLAGS) -fopenmp -o $@ $<

lifeCycle_transition.o lifeCycle_transition.mod: lifeCycle_transition.f90 lifeCycle_solveDP.o lifeCycle_algs.o share.f90
	$(FC) -c $(FFLAGS) -fopenmp -o $@ $<

lifeCycle_solveDP.o lifeCycle_solveDP.mod: lifeCycle_solveDP.f90 lifeCycle_algs.o share.o lifeCycle_vfuncs.o
	$(FC) -c $(FFLAGS) -fopenmp -o $@ $<

lifeCycle_algs.o lifeCycle_algs.mod: lifeCycle_algs.f90 lifeCycle_vfuncs.o share.o nrtype.o nrutil.o nr.o
	$(FC) -c $(FFLAGS) -fopenmp -o $@ $<

lifeCycle_vfuncs.o lifeCycle_vfuncs.mod: lifeCycle_vfuncs.f90 share.o
	$(FC) -c $(FFLAGS) -o $@ $<

share.o share.mod: shareIter.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%util.o %util.mod:  %util.f90 nrtype.o
	$(FC) -c $(FFLAGS) -o $@ $<

nr.o nr.mod: nr.f90 nrtype.o
	$(FC) -c $(FFLAGS) -o $@ $<

nrtype.o nrtype.mod: nrtype.f90
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod
