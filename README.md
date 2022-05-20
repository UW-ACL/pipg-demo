# PIPG: Proportional-Integral Projected Gradient Algorithm for Trajectory Optimization

- `./PIPG` contains the implementaion of PIPG for conic optimization as described in [Yu et al. 2022](https://doi.org/10.1016/j.automatica.2022.110359). The two examples in the paper are formulated in `ex1_problem_data.jl` and `ex2_problem_data.jl`. The simulations are carried out in the corresponding notebooks `ex1_problem_solve.jl` and `ex2_problem_solve.jl`.
- `PIPGeq_demo.ipynb` implements an earlier version of PIPG (called PIPGeq) described in [Yu et al. 2020](https://doi.org/10.1109/LCSYS.2020.3044977)
- `./exPIPG` contains the implementation of the _extrapolated_ PIPG for infeasibility detection. 
