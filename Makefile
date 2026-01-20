
MPIF90=mpifrtpx
OPTFLG=-Kfast
OMPFLG=-Kopenmp
F90=frtpx

all: euler_maruyama record_gen montecalro_w_record

euler_maruyama: src/double_well_SingleFil.f90
	${F90} ${OPTFLG} ${OMPFLG} src/double_well_SingleFil.f90 -o $@

record_gen: src/double_well_SingleFil_dictgen.f90
	${MPIF90} ${OPTFLG} ${OMPFLG} src/double_well_SingleFil_dictgen.f90 -o $@

montecalro_w_record: src/double_well_SingleFil_mod.f90
	${F90} ${OPTFLG} ${OMPFLG} src/double_well_SingleFil_mod.f90 -o $@

.PHONY: clean
clean:
	rm -f euler_maruyama record_gen montecalro_w_record *.o *.mod
	@echo "Cleaned up binary and module files."