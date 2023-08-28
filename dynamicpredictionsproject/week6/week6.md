# Week 6 Work
Arin Vansomphone

NSF Biostats Internship Summer 2023

## Brief Overview of Clinical Trial Program

### The Problem

![clinicaltrialimage.png](attachment:clinicaltrialimage.png)

* Clinical Trials have become increasingly more difficult to conduct (staffing, budget, and protocol concerns)
* More than 20% of clinical trials fail to recruit enough patients in time
* Biomedical companies (and ourselves!) rely on trial data to progress their treatment development
* How can we optimize multi-center, multi-state recruitment designs for global clinical trials?

### Underlying Mathematics: An Optimization Problem

* Recruit patients from $J$ countries (Anisimov, 2020)
* For a center in country $j$:

![Diagram%20of%20N%20Distribution.png](attachment:Diagram%20of%20N%20Distribution.png)

* Mean enrollment rate $m_j$, variance of enrollment rate $s^2_j$

**Minimize cost, constrainted on a probability of success, and solution must be integer - mixed integer non-linear programming**

### Example workflow


```julia
using ClinicalTrialOptm
using Distributions, HiGHS, SCS, Pajarito, SCIP, MathOptInterface, JuMP, Plots, StatsPlots
const MOI = MathOptInterface;
```


```julia
m = [1, 1.2, 1.4, 1.9] # mean enrollment rates per center
s² = [1.5, 1.7, 1.3, 1.1] # variances of enrollment rates
l = [0, 4, 2, 1] # lower bound of centers
u = [10, 24, 20, 15] # upper bound of center 
c₀ = [19000.0, 15000.0, 14000.0, 16000.0] # cost of initializing one center
c = [7000.0, 5000.0, 5000.0, 6000.0] # cost of running one center per unit of time 
q = [1000.0, 2000.0, 1500.0, 1600.0] # cost per one enrolled patient
d = [0.10, 0.09, 0.04, 0.07] # dropout rate
T₀ = fill(Uniform(0.0, 6.0), 4) # initialization period
Td = 12.0 # total time of clinical trial
ct = ClinicalTrial(m, s², l, u, c₀, c, q, d, T₀, Td)
```




    
    Global Clinical Trial:
    
    Optimal center assignment not computed.
    ┌─────────┬─────────┬────────┬────────────┬────────────────┬──────────────┬───────────────┬──────────────┬──────────────┬─────────┐
    │[1m Country [0m│[1m mean(λ) [0m│[1m var(λ) [0m│[1m init. cost [0m│[1m    maint. cost [0m│[1m enroll. cost [0m│[1m drop out rate [0m│[1m min. centers [0m│[1m max. centers [0m│[1m centers [0m│
    │[90m         [0m│[90m         [0m│[90m        [0m│[90m   $/center [0m│[90m $/center/month [0m│[90m    $/patient [0m│[90m               [0m│[90m              [0m│[90m              [0m│[90m         [0m│
    ├─────────┼─────────┼────────┼────────────┼────────────────┼──────────────┼───────────────┼──────────────┼──────────────┼─────────┤
    │       1 │    1.00 │   1.50 │      19000 │           7000 │         1000 │          0.10 │            0 │           10 │      NA │
    │       2 │    1.20 │   1.70 │      15000 │           5000 │         2000 │          0.09 │            4 │           24 │      NA │
    │       3 │    1.40 │   1.30 │      14000 │           5000 │         1500 │          0.04 │            2 │           20 │      NA │
    │       4 │    1.90 │   1.10 │      16000 │           6000 │         1600 │          0.07 │            1 │           15 │      NA │
    └─────────┴─────────┴────────┴────────────┴────────────────┴──────────────┴───────────────┴──────────────┴──────────────┴─────────┘
    
    





```julia
optdes!(ct, 400, ps = 0.95)
```

    Lower bound probability of success: 7.112034326831168e-28
    The optimal solution is not the lower bound of the centers.
    Upper bound probability of success: 0.9999940867534107
    The optimal solution is feasible.
      2.437051 seconds (6.64 M allocations: 356.935 MiB, 1.76% gc time, 98.88% compilation time: 25% of which was recompilation)
    solution_summary(model) = * Solver : SCIP
    
    * Status
      Result count       : 7
      Termination status : OPTIMAL
      Message from the solver:
      "SCIP_STATUS_OPTIMAL"
    
    * Candidate solution (result #1)
      Primal status      : FEASIBLE_POINT
      Dual status        : NO_SOLUTION
      Objective value    : 3.27739e+06
      Objective bound    : 3.27739e+06
      Relative gap       : 0.00000e+00
    
    * Work counters
      Solve time (sec)   : 1.78320e-02
      Simplex iterations : 89
      Node count         : 28
    
    termination_status(model) = MathOptInterface.OPTIMAL
    primal_status(model) = MathOptInterface.FEASIBLE_POINT
    objective_value(model) = 3.277387200000001e6





    
    Global Clinical Trial:
    
    Optimal center assignment calculated.
    An optimal solution has been found.
    Total duration (months): 12.0
    Target enrollment: 400
    Probability of success (based on normal approximation): 0.9550669693996681
    Probability of success (based on Poisson-Gamma model): 0.9636054293917478
    Expected cost ($): 3.2773872e6
    ┌─────────┬─────────┬────────┬────────────┬────────────────┬──────────────┬───────────────┬──────────────┬──────────────┬─────────┐
    │[1m Country [0m│[1m mean(λ) [0m│[1m var(λ) [0m│[1m init. cost [0m│[1m    maint. cost [0m│[1m enroll. cost [0m│[1m drop out rate [0m│[1m min. centers [0m│[1m max. centers [0m│[1m centers [0m│
    │[90m         [0m│[90m         [0m│[90m        [0m│[90m   $/center [0m│[90m $/center/month [0m│[90m    $/patient [0m│[90m               [0m│[90m              [0m│[90m              [0m│[90m         [0m│
    ├─────────┼─────────┼────────┼────────────┼────────────────┼──────────────┼───────────────┼──────────────┼──────────────┼─────────┤
    │       1 │    1.00 │   1.50 │      19000 │           7000 │         1000 │          0.10 │            0 │           10 │       0 │
    │       2 │    1.20 │   1.70 │      15000 │           5000 │         2000 │          0.09 │            4 │           24 │       5 │
    │       3 │    1.40 │   1.30 │      14000 │           5000 │         1500 │          0.04 │            2 │           20 │      20 │
    │       4 │    1.90 │   1.10 │      16000 │           6000 │         1600 │          0.07 │            1 │           15 │      14 │
    └─────────┴─────────┴────────┴────────────┴────────────────┴──────────────┴───────────────┴──────────────┴──────────────┴─────────┘
    
    




## Recent Progress

* Changed main solver from KNITRO (commercial) to SCIP (open-source)
* Created a documentation page on GitHub: https://hua-zhou.github.io/ClinicalTrialOptm.jl/dev/
* Initialized unit testing for package
* Working on binder demonstration page
* Plan on setting up code coverage website in future
