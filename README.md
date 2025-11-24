# LyapOpt vs MaxWeight — MATLAB Code

This repository contains the MATLAB code used in the paper:

**Finite‑Time Minimax Bounds and an Optimal Lyapunov Policy in Queueing Control**  
Yujie Liu, Vincent Y. F. Tan, Yunbei Xu  

# main files
- genDepartureSet.m # Generate the scheduling set \( \mathcal{D} \) and its facets.  
- genLambda.m # Generate random arrival rate vectors \( \lambda \) within the capacity region \( \Pi(\mathcal{D}) \). 
- hard_instance_b_change.m # Plot the \( \sqrt{b} \)-gap for hard instance used in Corollary 1.  
- hard_instance_t_change.m # Plot the \( \sqrt{t} \)-gap for hard instance used in Corollary 1.  
- run_single_lambda_plot.m # Plot total queue length/sum of squares of queue lengths for a single instance (lambda)
- run_single_lambda_before_arrival_plot.m # Plot total queue length/total queue lengths before arrivals for a single instance (lambda)
- run_lambda_performance.m # Batch performance evaluation over many instances (lambda)
- output_genDepartureSet_n=7.mat # Pre-generated D and faces for n=7
- output_genDepartureSet_n=8.mat # Pre-generated D and faces for n=8
