set pagination off
set confirm off
set disable-randomization off
set debuginfod enabled off
set auto-load python-scripts off
break correction_term_
run
set language c
printf "PARAMETER_READBACK xi_solv=%.17g eps=%.17g class=%c\n", *(double*) &__parametry_MOD_xi_solv, *(double*) &__parametry_MOD_eps, *(char*) &__parametry_MOD_solvclass
printf "UNIT_READBACK unitenergy=%.17g unitlength=%.17g kappaH=%.17g pH=%.17g zeta=%.17g\n", *(double*) &__lxlylz_MOD_unitenergy, *(double*) &__lxlylz_MOD_unitlength, *((double*) &__parametry_MOD_kappa_a), *((double*) &__parametry_MOD_p_a), *(double*) &__parametry_MOD_zeta
quit
