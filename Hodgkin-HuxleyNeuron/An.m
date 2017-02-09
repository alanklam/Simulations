function [y] = An(x)
y = (0.55 + 0.01*x)/(-exp(-5.5 - 0.1*x) + 1);