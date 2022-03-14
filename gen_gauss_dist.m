function [x, fx] = gen_gauss_dist(mu,sigma,interval,step)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%interval should be something like [-5,5]

x = interval(1):step:interval(2);
fx = (1/sigma*sqrt(2*pi))*exp(-(x-mu).^2/(2*sigma^2));

plot(x,fx);
grid on;

end

