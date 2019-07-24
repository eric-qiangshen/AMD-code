# AMD - Toolbox for computing separating input in set-membership identification of affine systems

This is the code of active model discrimination that computes the separating input for set-membership identification of affine systems. Please refer to the reference: 

Ding, Y., Harirchi F., Yong, S. Z., Jacobsen, E., and Ozay, N. (2018). Optimal input design for affine model discrimination with applications in intention-aware vehicles, in Proceedings of the 9th ACM/IEEE International Conference on Cyber-Physical Systems, PP. 297-307, arXiv:1702.01112.

## List of capabilities: 

 - Concantenating the affine dynamics over a finite horizon and a finite number of model/member pairs
 - Convert semi-infinite constraint to linear constraint by leveraging robust optimization
 - Implement double negation to find equivalence of the (non-convex) separability condition 
 - Solve a mixed integer linear program (MILP) to get the separating input for set-membership identification
   - Polyhedral constraints on input, controlled state, noise
 - Compare with previous work 
   - Harirchi F., Yong, S. Z., Jacobsen E., and Ozay, N. (2017). Active mode discrimination with applications to fraud  
     detection in smart buildings. IFAC-PapersOnLine, 50(1):9527-9534.

## Requirements
 - Recent version of Matlab
 - Optimizaiton modeling toolbox, Yalmip (https://yalmip.github.io/download/)
 - Solver: Gurobi for solving LP and supporting SOS-1 constraint (https://www.gurobi.com/)

