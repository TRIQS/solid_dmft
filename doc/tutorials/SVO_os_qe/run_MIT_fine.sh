#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
for J in 0.0 0.25 0.5 0.75 1.0
  do
  echo "Creating directory for J=$J eV"
  for U in 4.25 4.5 4.75 5 5.25 5.5 5.75 6 6.25 6.5 6.75 7 7.25
    do 
       echo "Creating directory for U=$U eV"
       mkdir -p "./J$J/U$U/"
       cp ./svo.h5  "./J$J/U$U/svo.h5"
       cat > "./J$J/U$U/dmft_config.toml" << EOF
[general]
seedname = "svo"
jobname = "out"
enforce_off_diag = true
block_threshold = 0.001

prec_mu = 0.001

h_int_type = "kanamori"
U = $U
J = $J

mu_initial_guess = 12.297745

beta = 40

n_iter_dmft = 8

dc_type = 1
dc = true
dc_dmft = false

calc_energies = false
sigma_mix = 1.0

h5_save_freq = 2

[solver]
type = "cthyb"
n_l = 35
length_cycle = 120
n_warmup_cycles = 8000
n_cycles_tot = 10e+6
imag_threshold = 1e-5
measure_G_l = true
perform_tail_fit = false

EOF
  cd "./J$J/U$U/"
  echo "I am in folder:"
  pwd
  echo "Running the DMFT cycle for parameters: U=$U eV, J=$J eV"
  
  mpirun -n 72 solid_dmft
  
  cd ..
  cd ..
 done
 done
