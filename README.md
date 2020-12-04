# MotionPlanning_2D_m

This package contains motion planning examples that runs the following four solvers: PSGCFS, RRT*-CFS, CFS, CHOMP.

To run these examples (main files), users will need to install supporting packages.
1. Run `install_mpt3.m`.
2. Install other packeges suggested by MATLAB when running the examples.  

There are four examples in this package.
- `extendedCFS_original.m` is a simple implementation of the Convex Feasible Set (CFS) algorithm [1] in 2D planning.
- `main.m` implements CHOMP, CFS, and projected stochastic gradient Convex Feasible Set (PSGCFS) and compares the 2D planning performances of these methods.
- `main_MPC.m` implements CFS and PSGCFS in a MPC framework. The defult setting shows that PSGCFS can avoid local optima from which CFS cannot escape.
- `main_RRT_CFS_2.m` implements RRT*-CFS [2] for 2D planning. RRT*-CFS can solve narrow passage problems that are generally difficult for gradient-based solvers. 

![GitHub Logo](/RRT-CFS.jpg)

[1] Liu, Changliu, Chung-Yen Lin, and Masayoshi Tomizuka. "The convex feasible set algorithm for real time optimization in motion planning." SIAM Journal on Control and optimization 56.4 (2018): 2712-2733.

[2] https://jessicaleu24.github.io/ACC2021.html
