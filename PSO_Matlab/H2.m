function [value, normalized_infeasability] = H2(position)

x = position(1);
y = position(2);

%d = sqrt((x-8.6998)^2+(y-6.7665)^2);

value = 0.5 - ((sin(sqrt(x^2 + y^2)))^2 - 0.5)/(1 + 0.001*(x^2 + y^2))^2;
value =  - value; %to invert the problem to search for maxima 
%                     %with a minimization algorithm
normalized_infeasability = 0;