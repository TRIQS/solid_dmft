
[dft]: DFT related inputs
-------------------------

List of parameters that relate to the DFT calculation, useful mostly when doing CSC.






.. admonition:: dft_code 
 	:class: intag  
 
            **type=** string

            Choose the DFT code interface, for now Quantum Espresso and Vasp are available.

            Possible values:

            * 'vasp'
            * 'qe'

.. admonition:: n_cores 
 	:class: intag  
 
            **type=** int

            number of cores for the DFT code (VASP)

.. admonition:: n_iter 
 	:class: intag  
 
            **type=** int;  **optional**;  **default=**  6

            number of dft iterations per cycle

.. admonition:: dft_exec 
 	:class: intag  
 
            **type=** string;  **default=**  'vasp_std'

            command for the DFT / VASP executable

.. admonition:: store_eigenvals 
 	:class: intag  
 
            **type=** bool;  **optional**;  **default=**  False

            stores the dft eigenvals from LOCPROJ (projector_type=plo) or
            wannier90.eig (projector_type=w90) file in h5 archive

.. admonition:: mpi_env 
 	:class: intag  
 
            **type=** string;  **default=**  'local'

            selection for mpi env for DFT / VASP in default this will only call VASP as mpirun -np n_cores_dft dft_exec

.. admonition:: projector_type 
 	:class: intag  
 
            **type=** string;  **optional**;  **default=**  'plo'

            plo: uses VASP's PLO formalism, requires LOCPROJ in the INCAR
            w90: uses Wannier90, requires LWANNIER90 in the INCAR

.. admonition:: w90_exec 
 	:class: intag  
 
            **type=** string;  **default=** 'wannier90.x'

            the command to start a single-core wannier run

.. admonition:: w90_tolerance 
 	:class: intag  
 
            **type=** float;  **default=** 1e-6

            threshold for mapping of shells and checks of the Hamiltonian
