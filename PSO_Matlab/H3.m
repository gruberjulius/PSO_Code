function [value, normalized_infeasability] = H3(position)

n = 16; %dimensionality of problem
x = zeros(n,1);
t = 0.05;
s = 0.2;
c = 0.15;

for i=1:n,
x(i) = position(i);
end
value = 0;

for i=1:n,
    z_i = floor(abs(x(i)/s)+0.49999)*sign(x(i))*s;
    if rem(i,4)==1,
        d_i = 1;
    elseif rem(i,4)==2,
        d_i = 1000;
    elseif rem(i,4)==3,
        d_i = 10;
    else
        d_i = 100;
    end
    if abs(x(i)-z_i)<t,        
        value = value + (t*sign(z_i) + z_i)^2 * c * d_i;
    else 
        value = value + d_i*x(i)^2;
    end
end

normalized_infeasability = 0;
