[advanced]: Advanced inputs
------

Advanced parameters, do not modify default value unless you know what you are doing

List of possible entries:

dc_factor; dc_fixed_value; dc_fixed_occ; dc_U; dc_J; map_solver_struct; mapped_solver_struct_degeneracies; 


.. admonition:: dc_factor: 
 
            **type=** float;  **optional**;  **default=**  'none' (corresponds to 1)

            If given, scales the dc energy by multiplying with this factor, usually < 1

.. admonition:: dc_fixed_value: 
 
            **type=** float;  **optional**;  **default=**  'none'

            If given, it sets the DC (energy/imp) to this fixed value. Overwrites EVERY other DC configuration parameter if DC is turned on

.. admonition:: dc_fixed_occ: 
 
            **type=** list of float;  **optional**;  **default=**  'none'

            If given, the occupation for the DC for each impurity is set to the provided value.
            Still uses the same kind of DC!

.. admonition:: dc_U: 
 
            **type=** float or comma seperated list of floats;  **optional**;  **default=**  general_params['U']

            U values for DC determination if only one value is given, the same U is assumed for all impurities

.. admonition:: dc_J: 
 
            **type=** float or comma seperated list of floats;  **optional**;  **default=**  general_params['J']

            J values for DC determination if only one value is given, the same J is assumed for all impurities

.. admonition:: map_solver_struct: 
 
            **type=** dict;  **optional**;  **default=** no additional mapping

            Additional manual mapping of the solver block structure, applied
            after the block structure finder to all impurities.

.. admonition:: mapped_solver_struct_degeneracies: 
 
            **type=** list;  **optional**;  **default=** none

