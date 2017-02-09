function [I] = Ion(n,m,h,V,gl,gk,gna,Vl,Vk,Vna)
I = gl*(V+Vl) + gk*(n^4)*(V + Vk) + gna*(m^3)*h*(V + Vna);
end

