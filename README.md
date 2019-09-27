# NoisyNetworkReconstruction_SDDV
Simulation of noisy coupled networks and the subsequent reconstruction; Simulations and Code

reconstruction_network_simulation: contains matlab functions for simulating time series and reconstruction algorithms (method *2p_alg*: approximate activity f and time derivative cdot at t<tau<(t+delta t) with f=0.5*(f(t)+f(t+delta t) and xdot = 1/delta t * f(t+delta t) - f(t)
*3p_alg*: approximate time derivative at time t with f(t+delta t) - f(t-delta t), use activity at time t, 
*nextstep_alg*: reconstruct with f(x+ delta t) = A* f(x)

proof of concept: 
*strawman_wo_noise* simulate timeseries (one relaxation from random initial condition) without noise  and with noise, compare reconstructability
*precisiontest* for two different noise levels, compare reconstruction with 2p_alg for full precision, rounding after three digits and addition of noise with std (1/6 10^-3; then 99.8% of all values fall in 0.5 10^-3 ~ rounding)
*linear sigma* simulate tseries with different levels of noise, save SNR (based on expectation value) and reconstruction AUC for 2p and nextstep algorithm

optimal regime:
*vary sigma cubic* simulate tseries with cubic decay with different levels of noise, save SNR (based on expectation value) and reconstruction AUC for 2p, 3p and nextstep algorithm
*timescales* nrepeat above experiment with different timescales (decrease coarse scale, increase fine scale)


