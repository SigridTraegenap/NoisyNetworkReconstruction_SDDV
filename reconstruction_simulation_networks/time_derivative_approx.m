function [taus,x_tau,derivative]=time_derivative_approx(activity, T, delta_t)
%calculate  approximation of time derivatives with point in future approximation
%xdot = 1/h * (x(t+h) - x(t))
%is only valid at point tau
%estimate activity also at that point tau
time_points = (1:T/delta_t)*delta_t;
taus = (time_points(2:end)-time_points(1:end-1))/2;
derivative = (activity(:,2:end) - activity(:,1:end-1))/delta_t;
x_tau = (activity(:,2:end) + activity(:,1:end-1))/2;
end