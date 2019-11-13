# Main Code
This folder contains the main code that runs the algorithm
# B_Updated.m
This function applies a low rank Quasi newton approximation of the Hessian each step of the local transport function
# Global_optimal_transport.m
This function is the main, globally iterative code for the Weighted Adaptive Optimal Transport algorithm
# Local_transport_function.m
This function is run during each local step of the globally iterative transport function. It is called within a for loop in the global code and completes when the gradient is converged to zero between two nearby distributions of sampled points. 
# gh.m
This is the gradient function that depends on several functions within the Gradient_Dependencies folder. It calculates a standard gradient descent on two sets of points. 
