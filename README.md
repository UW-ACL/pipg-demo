<center> <h1>PIPG: Proportional-Integral Projected Gradient Algorithm for Trajectory Optimization 1</h1> </center>

<center> [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) </center>

- `./PIPG` contains the implementaion of PIPG for conic optimization as described in [Yu et al. 2022](https://doi.org/10.1016/j.automatica.2022.110359). The two examples in the paper are formulated in `ex1_problem_data.jl` and `ex2_problem_data.jl`. The simulations are carried out in the corresponding notebooks `ex1_problem_solve.jl` and `ex2_problem_solve.jl`.
- `PIPGeq_demo.ipynb` implements an earlier version of PIPG (called PIPGeq) described in [Yu et al. 2020](https://doi.org/10.1109/LCSYS.2020.3044977)
- `./exPIPG` contains the implementation of the _extrapolated_ PIPG for infeasibility detection. 

## Citations 

Conic optimization solver PIPGeq specialized to affine equality constraint and convex set constraints which allow efficient projections.
```
@article{pipgeq,
  author={Yu, Yue and Elango, Purnanand and Açıkmeşe, Behçet},
  journal={IEEE Control Systems Letters}, 
  title={Proportional-Integral Projected Gradient Method for Model Predictive Control}, 
  year={2021},
  volume={5},
  number={6},
  pages={2174-2179},
  doi={10.1109/LCSYS.2020.3044977}}
```

Conic optimization solver PIPG which can handle general conic constraints and convex set constraints which allow efficient projections.
```
@article{pipg,
  author = {Yue Yu and Purnanand Elango and Ufuk Topcu and Behçet Açıkmeşe},
  title = {Proportional–integral projected gradient method for conic optimization},
  journal = {Automatica},
  volume = {142},
  pages = {110359},
  year = {2022},
  issn = {0005-1098},
  doi = {https://doi.org/10.1016/j.automatica.2022.110359},
  url = {https://www.sciencedirect.com/science/article/pii/S0005109822002096}
}
```

Infeasibility detection using PIPG
```
@misc{pipg-infeas,
  author = {Yu, Yue and Topcu, Ufuk},    
  title = {Proportional-Integral Projected Gradient Method for Infeasibility Detection in Conic Optimization},
  publisher = {arXiv},  
  year = {2021},
  doi = {10.48550/ARXIV.2109.02756},  
  url = {https://arxiv.org/abs/2109.02756}
}

```

Extrapolated PIPG
```
@misc{expipg,
  author = {Yu, Yue and Elango, Purnanand and Açıkmeşe, Behçet and Topcu, Ufuk},  
  title = {Extrapolated Proportional-Integral Projected Gradient Method for Conic Optimization},
  publisher = {arXiv},
  year = {2022},
  doi = {10.48550/ARXIV.2203.04188},
  url = {https://arxiv.org/abs/2203.04188}  
}
```

Customized PIPG for real-time powered-descent guidance
```
@inproceedings{pipg-pdg,
  author={Elango, Purnanand and Kamath, Abhinav G. and Yu, Yue and Carson III, John M. and Mesbahi, Mehran and Açıkmeşe, Behçet}, 
  title={A Customized First-Order Solver for Real-Time Powered-Descent Guidance}, 
  booktitle={AIAA SciTech 2022 Forum}, 
  publisher={American Institute of Aeronautics and Astronautics}, 
  year={2022}, 
  month={Jan},  
  doi={10.2514/6.2022-0951}, 
  url={https://arc.aiaa.org/doi/abs/10.2514/6.2022-0951}
}

```
