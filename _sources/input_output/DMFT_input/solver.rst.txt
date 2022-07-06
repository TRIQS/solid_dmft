[solver]: solver specific parameters
------------------------------------

Here are the parameters that are uniquely dependent on the solver chosen. Below a list of the supported solvers:








.. admonition:: store_solver 
 	:class: intag  
 
            **type=** bool;  **optional** default= False

            store the whole solver object under DMFT_input in h5 archive

cthyb parameters
================

.. admonition:: length_cycle 
 	:class: intag  
 
            **type=** int

            length of each cycle; number of sweeps before measurement is taken

.. admonition:: n_warmup_cycles 
 	:class: intag  
 
            **type=** int

            number of warmup cycles before real measurement sets in

.. admonition:: n_cycles_tot 
 	:class: intag  
 
            **type=** int

            total number of sweeps

.. admonition:: measure_G_l 
 	:class: intag  
 
            **type=** bool

            measure Legendre Greens function

.. admonition:: measure_G_tau 
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=** True

            should the solver measure G(tau)?

.. admonition:: measure_G_iw 
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=** False

            should the solver measure G(iw)?

.. admonition:: measure_density_matrix 
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=** False

            measures the impurity density matrix and sets also
            use_norm_as_weight to true

.. admonition:: measure_pert_order 
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=** False

            measure perturbation order histograms: triqs.github.io/cthyb/latest/guide/perturbation_order_notebook.html

            The result is stored in the h5 archive under 'DMFT_results' at every iteration
            in the subgroups 'pert_order_imp_X' and 'pert_order_total_imp_X'

.. admonition:: max_time 
 	:class: intag  
 
            **type=** int;  **optional**;  **default=** -1

            maximum amount the solver is allowed to spend in each iteration

.. admonition:: imag_threshold 
 	:class: intag  
 
            **type=** float;  **optional**;  **default=**  10e-15

            threshold for imag part of G0_tau. be warned if symmetries are off in projection scheme imag parts can occur in G0_tau

.. admonition:: off_diag_threshold 
 	:class: intag  
 
            **type=** float;  **optional**

            threshold for off-diag elements in Hloc0

.. admonition:: delta_interface 
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=** False

            use new delta interface in cthyb instead of input G0

.. admonition:: move_double 
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=** True

            double moves in solver

.. admonition:: perform_tail_fit 
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=** False

            tail fitting if legendre is off?

.. admonition:: fit_max_moment 
 	:class: intag  
 
            **type=** int;  **optional**

            max moment to be fitted

.. admonition:: fit_min_n 
 	:class: intag  
 
            **type=** int;  **optional**

            number of start matsubara frequency to start with

.. admonition:: fit_max_n 
 	:class: intag  
 
            **type=** int;  **optional**

            number of highest matsubara frequency to fit

.. admonition:: fit_min_w 
 	:class: intag  
 
            **type=** float;  **optional**

            start matsubara frequency to start with

.. admonition:: fit_max_w 
 	:class: intag  
 
            **type=** float;  **optional**

            highest matsubara frequency to fit

.. admonition:: random_seed 
 	:class: intag  
 
            **type=** int;  **optional** default by triqs

            if specified the int will be used for random seeds! Careful, this will give the same random
            numbers on all mpi ranks

.. admonition:: legendre_fit 
 	:class: intag  
 
            **type=** bool;  **optional** default= False

            filter noise of G(tau) with G_l, cutoff is taken from n_l

ftps parameters
===============

.. admonition:: n_bath 
 	:class: intag  
 
            **type=** int

            number of bath sites

.. admonition:: bath_fit 
 	:class: intag  
 
            **type=** bool;  **default=** False

            DiscretizeBath vs BathFitter

.. admonition:: refine_factor 
 	:class: intag  
 
            **type=** int;  **optional**;  **default=** 1

            rerun ftps cycle with increased accuracy

.. admonition:: ph_symm 
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=** False

            particle-hole symmetric problem

.. admonition:: calc_me 
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=** True

            calculate only symmetry-inequivalent spins/orbitals, symmetrized afterwards

.. admonition:: enforce_gap 
 	:class: intag  
 
            **type=** list of floats;  **optional**;  **default=** 'none'

            enforce gap in DiscretizeBath between interval

.. admonition:: ignore_weight 
 	:class: intag  
 
            **type=** float;  **optional**;  **default=** 0.0

            ignore weight of peaks for bath fitter

.. admonition:: dt 
 	:class: intag  
 
            **type=** float

            time step

.. admonition:: state_storage 
 	:class: intag  
 
            **type=** string;  **default=**  './'

            location of large MPS states

.. admonition:: path_to_gs 
 	:class: intag  
 
            **type=** string;  **default=**  'none'

            location of GS if already present. Use 'postprocess' to skip solver and go directly to post-processing
            of previously terminated time-evolved state

.. admonition:: sweeps 
 	:class: intag  
 
            **type=** int;  **optional**;  **default=**  10

            Number of DMRG sweeps

.. admonition:: maxmI 
 	:class: intag  
 
            **type=** int;  **optional**;  **default=**  100

            maximal imp-imp bond dimensions

.. admonition:: maxmIB 
 	:class: intag  
 
            **type=** int;  **optional**;  **default=**  100

            maximal imp-bath bond dimensions

.. admonition:: maxmB 
 	:class: intag  
 
            **type=** int;  **optional**;  **default=**  100

            maximal bath-bath bond dimensions

.. admonition:: tw 
 	:class: intag  
 
            **type=** float, default 1E-9

            truncated weight for every link

.. admonition:: dmrg_maxmI 
 	:class: intag  
 
            **type=** int;  **optional**;  **default=**  100

            maximal imp-imp bond dimensions

.. admonition:: dmrg_maxmIB 
 	:class: intag  
 
            **type=** int;  **optional**;  **default=**  100

            maximal imp-bath bond dimensions

.. admonition:: dmrg_maxmB 
 	:class: intag  
 
            **type=** int;  **optional**;  **default=**  100

            maximal bath-bath bond dimensions

.. admonition:: dmrg_tw 
 	:class: intag  
 
            **type=** float, default 1E-9

            truncated weight for every link

ctseg parameters
================

.. admonition:: measure_hist 
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=** False

               measure perturbation_order histograms

.. admonition:: improved_estimator  
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=** False

              measure improved estimators
              Sigma_iw will automatically be calculated via
              http://dx.doi.org/10.1103/PhysRevB.85.205106
