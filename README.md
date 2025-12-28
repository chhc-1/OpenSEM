# Synthetic_eddy_cpp

WIP library for testing various Synthetic Eddy methods for turbulent inflow. Synthetic eddy methods currently in package (and may not be complete):

- Original Synthetic Eddy method (oSEM) (Jarrin et al. 2006) [1]

- Improved Synthetic Eddy method (ISEM) (Xiong et al. 2024) [2]

- Multi-region Synthetic eddy method (MRSEM) (Pamiès et al. 2009) [3]

- Divergence-free Synthetic eddy method (DFSEM) (Poletto 2013) [4]

Note that the DFSEM implemented has only been tested with an isotropic Reynold's Stress Tensor and only works with a very limited range of Reynold's Stress Tensors. (Corresponds to DFSEMiso from [4])

### References for Synthetic eddy methods

[1] N. Jarrin, S. Benhamadouche, D. Laurence, and R. Prosser. A synthetic-eddy-method for generating inflow conditions for large-eddy simulations. International Journal of Heat and Fluid Flow, vol. 27, no. 4, pp. 585–593, 2006.

[2] D. Xiong, Y. Yang, and Y. Wang. An improved synthetic eddy method for generating inlet turbulent boundary layers. Aerospace, vol. 9, no. 1, pp. 37–37, 2022.

[3] M. Pamiès, P.-É. Weiss, É. Garnier, S. Deck, and P. Sagaut. Generation of synthetic turbulent inflow data for large eddy simulation of spatially evolving wall-bounded flows. Physics of Fluids, vol. 21, no. 4, 2009.

[4] R. Poletto. Divergence free development of the syntehtic eddy method in order to improve synthetic turbulence for embedded LES simulations. PhD thesis, University of Manchester, 2013.