function [derivative]=time_derivative_IMapprox(activity, delta_t)
%calculate  approximation of time derivatives with two point approximation
%xdot = 1/2h * (x(t+h) - x(t-h))
derivative = (activity(:,3:end) - activity(:,1:end-2))/(2*delta_t);
end