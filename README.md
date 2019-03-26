# NoisyNetworkReconstruction_SDDV
Simulation of noisy coupled networks and the subsequent reconstruction; Simulations and Code

reconstruction_network_simulation: contains simple matlab function for simulating time series and reconstruction algorithms (method 1: approximate time derivative at tau f(x+h) - f(x), 2: approximate at t with f(x+h) - f(x-h), 23: reconstruct with f(x+ delta) = A* f(x)

proof of concept: simulate time series, 
simulation1: contains all simulated tseries and summary statistics for different algorithms (summary_...)

different timescales: containes cubic decay tseries simualted with dt=10-5 or dt=10-3

