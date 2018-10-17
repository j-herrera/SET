SET
===

Overview
--------
SET (Shock Estimation Tool) is the implementation of an approximate method for shock structure estimation in supersonic and hypersonic vehicles.

The procedure applies a combination of supersonic flow theories over surface streamlines, generating local shocks. These are then meshed together to obtain one or multiple global shocks, taking into account multiple origins through stagnation point clustering. Detached shocks are
addressed starting the procedure from an engineering correlation around the stagnation point and correcting the expansion through a correction factor.

More information on the methodology and some examples of application can be found in:

[Shock estimation in supersonic vehicles](https://suprimo.lib.strath.ac.uk:443/SUVU01:LSCOP_SU:SUDIGI28628)

[Shock-conforming mesh generation for aerodynamic analyses at supersonic regimes](https://www.sciencedirect.com/science/article/abs/pii/S0045793017303250)

Building
--------
SET is a python module that wraps a C++ shared library. The library is under cgal_wrapper and has to be compiled before use. CGAL, Python 3, Boost and Eigen3 are dependencies.

The build steps are as follows:

1. Execute 'cmake .' inside the cgal_wrapper folder to generate a makefile
2. Execute 'make' to compile the library into wrapper.so
3. Move wrapper.so into the SET folder that contains the python code