
Iterations
----------

List of the main outputs for solid_dmft for every iteration.

.. warning::

  According to the symmetries found by the solver, the resulting indexing of the triqs.Gf objects might vary.
  In order to retrieve the indices call the Gf.indices method.


Legend:

* iiter = iteration number: range(0, n_dmft_iter)
* ish = shell number: range(0, n_shells)
* icrsh = correlated shell number: range(0, n_corr_shells)
* iineq = inequivalent correlated shell number: range(0, n_inequiv_shells)
* iorb = orbital number: range(0, n_orbitals)
* sp = spin label
* ikpt = k-point label, the order is the same as given in the wannier90 input: range(0, n_kpt)
* iband = band label before downfolding, n_bands = number of bands included in the disentanglement window during the wannierization: range(0, n_bands)


[observables]
=============

.. admonition:: chemical_potential_pre: 
  :class: intag
            **type=** float;

            Chemical potential before the solver iteration.

.. admonition:: chemical_potential_post: 
  :class: intag
            **type=** float;

            Chemical potential after the solver iteration.

.. admonition:: DC_energ: 
  :class: intag
            **type=** arr(float);

            **indices=** [iorb]

            Double counting correction.

.. admonition:: DC_pot: 
  :class: intag
 
            **type=** arr(float);

            **indices=** [iiter]

            Double counting potential.**what exactly is the indexing here?**

.. admonition:: Delta_time_{iimp}: 
  :class: intag
 
            **type=** triqs.gf.block_gf.BlockGf


            Imaginary time hybridization function.

.. admonition:: G0_freq_{iimp}: 
  :class: intag
 
            **type=** triqs.gf.block_gf.BlockGf


            Imaginary frequency Weiss field.

.. admonition:: G0_time_orig_{iimp}: 
  :class: intag
 
            **type=** triqs.gf.block_gf.BlockGf


            ??

.. admonition:: G_imp_freq_{iimp}: 
  :class: intag
 
            **type=** triqs.gf.block_gf.BlockGf


            Imaginary frequency impurity green function.

.. admonition:: G_imp_l_{iimp}: 
  :class: intag
 
            **type=** triqs.gf.block_gf.BlockGf


            Legendre representation of the impurity green function.

.. admonition:: G_imp_time_{iimp}: 
  :class: intag
 
            **type=** triqs.gf.block_gf.BlockGf


            Imaginary time representation of the impurity green function.

.. admonition:: Sigma_freq_{iimp}: 
  :class: intag
 
            **type=** triqs.gf.block_gf.BlockGf


            Imaginary frequency self-energy obtained from the Dyson equation.

.. admonition:: deltaN: 
  :class: intag
 
            **type=** dict(arr(float))
            
            **indices=** [ispin][ikpt][iband, iband]


            Correction to the DFT occupation of a particular band: 

.. admonition:: deltaN_trace: 
  :class: intag
 
            **type=** dict
            
            **indices=** [ispin]


            Total sum of the charge correction for an impurity.

.. admonition:: dens_mat_pre: 
  :class: intag
 
            **type=** arr(dict) 

            **indices=** [iimp][*same as block structure Gf*]

            Density matrix before the solver iteration.

.. admonition:: dens_mat_post: 
  :class: intag
 
            **type=** arr(dict) 

            **indices=** [ispin][iimp]

            Density matrix after the solver iteration.

