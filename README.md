# Analytical calculation of thermoelastic stress
Calculating material response due to rapid heating from laser pulse. The analytical solution is obtained by resolving the laser pulse as a Fourier series and calculating the the electronic and lattice temperature using the two-temperature model. Stress, temperature, and displacement field are outputted as a function of the desired laser fluence and pulse duration.

TECHNICAL DETAILS:

        VERSIONS:

                - MATLAB (Codes were tested with MATLAB R2024A)
TASKS:

        1) analytical_thermo_elastic_stress

                This code is to obtain the stress, temperature, and displacement field by solving the two-temperature model coupled to an equation of motion.
		Here, the user can input the material properties, and select the range (or single value) of the desired laser fluence and pulse duration.

        2) compute_coefficients
                
                This code reads the pulse duration value inputted in **analytical_thermo_elastic_stress** and outputs the fourier coefficients back to 
		**analytical_thermo_elastic_stress** to obtain the stress, temperature, and displacement field.
                
		Please see [[CITE]](https://arxiv.org/abs/2412.04762) for further details.
