# Proportional-Integral Projected Gradient Method (PIPG) for Trajectory Optimization

- `./pipg` contains the implementaion of PIPG for conic optimization as described in https://arxiv.org/abs/2108.10260. The two examples in the paper are formulated in `ex1_problem_data.jl` and `ex2_problem_data.jl`. The simulations are carried out in the corresponding notebooks `ex1_problem_solve.jl` and `ex2_problem_solve.jl`.
- `pipgeq_demo.jl` implements an earlier version of PIPG (called PIPGeq) described in https://doi.org/10.1109/LCSYS.2020.3044977
