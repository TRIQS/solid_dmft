
Observables/convergence_obs
---------------------------

List of the single-particle observables obtained in a single DMFT iteration


Legend:

* iiter = iteration number: range(0, n_dmft_iter)
* iimp = impurity number: range(0, n_imp)
* iorb = orbital number: range(0, n_orbitals)
* ispin = spin label, 'up' or 'down' in collinear calculations


[observables]
=============

.. admonition:: iteration:
  :class: intag

            **type=** arr(int);

            **indices=** [iiter]

            Number of the iteration.

.. admonition:: mu:
  :class: intag

            **type=** arr(float);

            **indices=** [iiter]

            Chemical potential fed to the solver at the present iteration (pre-dichotomy adjustment).

.. admonition:: orb_gb2:
  :class: intag

            **type=** arr(dict)

            **indices=** [iimp][ispin][iiter, iorb]

            Orbital resolved G(beta/2), proxy for projected density of states at the Fermi level. Low value of orb_gb2 correlate with the presence of a gap.

.. admonition:: imp_gb2:
  :class: intag

            **type=** arr(dict)

            **indices=** [iimp][ispin][iiter]

            Site G(beta/2), proxy for total density of states at the Fermi level. Low values correlate with the presence of a gap.

.. admonition:: orb_Z:
  :class: intag

            **type=** arr(dict)

            **indices=** [iimp][ispin][iiter, iorb]

            Orbital resolved quasiparticle weight (eff_mass/renormalized_mass). As obtained by linearizing the self-energy around :math:`\omega = 0`

            .. math::

              Z = \bigg( 1- \frac{\partial Re[\Sigma]}{\partial \omega} \bigg|_{\omega \rightarrow 0} \bigg)^{-1} \\


.. admonition:: orb_occ:
  :class: intag

            **type=** arr(dict)

            **indices=** [iimp][ispin][iiter, iorb]

            Orbital resolved mean site occupation.

.. admonition:: imp_occ:
  :class: intag

            **type=** arr(dict)

            **indices=** [iimp][ispin][iiter]

            Total mean site occupation.


.. admonition:: E_tot:
  :class: intag

            **type=** arr(float)

            **indices=** [iiter]

            Total energy, computed as:

            .. math::

              E_{tot} = E_{DFT} + E_{corr} + E_{int} -E_{DC}


.. admonition:: E_dft:
  :class: intag

            **type=** arr(float)

            **indices=** [iiter]

            :math:`E_{DFT}` in the total energy expression. System energy as computed by the DFT code at every csc iteration.



.. admonition:: E_bandcorr:
  :class: intag

            **type=** arr(float)

            **indices=** [iiter]

            :math:`E_{corr}` in the total energy expression. DMFT correction to the kinetic energy.

.. admonition:: E_corr_en:
  :class: intag

            **type=** arr(float)

            **indices=** [iiter]

            Sum of the E_DC and E_int_imp terms.

.. admonition:: E_int_imp:
  :class: intag

            **type=** arr(float)

            **indices=** [iiter]

            :math:`E_{int}` in the total energy expression. Energy contribution from the electronic interactions within the single impurity.


.. admonition:: E_DC:
  :class: intag

            **type=** arr(float)

            **indices=** [iiter]

            :math:`E_{DC}` in the total energy expression. Double counting energy contribution.




[convergence_obs]
=================

.. admonition:: iteration:
  :class: intag

            **type=** arr(int);

            **indices=** [iiter]

            Number of the iteration.

.. admonition:: d_mu:
  :class: intag

            **type=** arr(float)

            **indices=** [iiter]

            Chemical potential stepwise difference.


.. admonition:: d_orb_occ:
  :class: intag

            **type=** arr(dict)

            **indices=** [iimp][ispin][iiter,iorb]

            Orbital occupation stepwise difference.

.. admonition:: d_imp_occ:
  :class: intag

            **type=** arr(dict)

            **indices=** [iimp][ispin][iiter]

            Impurity occupation stepwise difference.

.. admonition:: d_Gimp:
  :class: intag

            **type=** arr(float)

            **indices=** [iiter]

            DMFT self-consistency condition | Gloc - Gimp | difference of current iteration.

.. admonition:: d_G0:
  :class: intag

            **type=** arr(float)

            **indices=** [iiter]

            Weiss field stepwise difference.

.. admonition:: d_Sigma:
  :class: intag

            **type=** arr(float)

            **indices=** [iiter]

            Impurity self-energy stepwise difference.


.. admonition:: d_Etot:
  :class: intag

            **type=** arr(float)

            **indices=** [iiter]

            Total energy stepwise difference.


