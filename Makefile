# 
# Purpose: Build new subcolumn generator
# Author: Benjamin R. Hillman
#

F90 = gfortran

subcol: gen_subcol.F90 gen_subcol_cld.F90 gen_subcol_var.F90 \
		adjust_precip.F90 \
        init_random_seed.F90 random.o quick_sort.f90 scops.f90 prec_scops.f90
	f2py --f90flags='-g -fbounds-check' -c -m $@ $^

# compile object files
%.o: %.F90
	$(F90) -fPIC -c $<

%.o: %.f90
	$(F90) -fPIC -c $<

%.o: %.f
	$(F90) -fPIC -c $<

# clean up
clean:
	rm -f *.o *.pyf *.mod *.so
