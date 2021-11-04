function [value, normalized_infeasability] = H1(position)

x = position(1);
y = position(2);

d = sqrt((x-8.6998)^2+(y-6.7665)^2);

value = ((sin(x-y/8))^2 + (sin(y+x/8))^2)/(d+1);
value = -value; %to invert the problem to search for maxima
normalized_infeasability = 0;

%problem has a maximum at (8.6698,6.7665) of 2.0
