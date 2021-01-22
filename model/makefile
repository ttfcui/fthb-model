
FC = gfortran
FFLAGS = -g -w -O3 -ffree-line-length-none -fcheck=bounds -fmax-errors=1
MODFLAGS = -fmax-stack-var-size=65535 -freal-4-real-8

SRCS = $(patsubst %.f90, %.o, $(wildcard lifeCycle*.f90) share.f90 nrutil.f90 nrtype.f90 cdflib.f90)


lifecycle.out: $(SRCS)
	$(FC) $(FFLAGS) -fopenmp -o $@ $^

lifeCycle_TC10172016.o: lifeCycle_TC10172016.f90 lifeCycle_models.o
	$(FC) -c $(FFLAGS) $(MODFLAGS) -fopenmp -o $@ $<

lifeCycle_models.o lifeCycle_models.mod: lifeCycle_models.f90 lifeCycle_simulate.o lifeCycle_transition.o
	$(FC) -c $(FFLAGS) $(MODFLAGS) -fopenmp -o $@ $<

lifeCycle_simulate.o lifeCycle_simulate.mod: lifeCycle_simulate.f90 lifeCycle_solveDP.o lifeCycle_algs.o lifeCycle_calibrate.o
	$(FC) -c $(FFLAGS) $(MODFLAGS) -fopenmp -o $@ $<

lifeCycle_transition.o lifeCycle_transition.mod: lifeCycle_transition.f90 lifeCycle_solveDP.o lifeCycle_algs.o lifeCycle_calibrate.o
	$(FC) -c $(FFLAGS) $(MODFLAGS) -fopenmp -o $@ $<

lifeCycle_solveDP.o lifeCycle_solveDP.mod: lifeCycle_solveDP.f90 lifeCycle_algs.o share.o lifeCycle_vfuncs.o
	$(FC) -c $(FFLAGS) $(MODFLAGS) -fopenmp -o $@ $<

lifeCycle_vfuncs.o lifeCycle_vfuncs.mod: lifeCycle_vfuncs.f90 lifeCycle_algs.o share.o
	$(FC) -c $(FFLAGS) $(MODFLAGS) -o $@ $<

lifeCycle_calibrate.o lifeCycle_calibrate.mod: lifeCycle_calibrate.f90 cdflib.o share.o lifeCycle_algs.o
	$(FC) -c $(FFLAGS) $(MODFLAGS) -o $@ $<

lifeCycle_algs.o lifeCycle_algs.mod: lifeCycle_algs.f90 share.o nrtype.o nrutil.o
	$(FC) -c $(FFLAGS) $(MODFLAGS) -fopenmp -o $@ $<

share.o share.mod: shareIter.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%util.o %util.mod:  %util.f90 nrtype.o
	$(FC) -c $(FFLAGS) -o $@ $<

cdflib.o cdflib.mod: cdflib.f90
	$(FC) -c $(FFLAGS) -o $@ $<

nrtype.o nrtype.mod: nrtype.f90
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	rm -f *.o *.mod
