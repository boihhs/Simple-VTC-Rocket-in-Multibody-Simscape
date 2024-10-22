m = 4; % mass in Kg
horizon = 100;
h = 1; % height in m
g = 9.81; % gravity constant
Radius = 95; % Radius in m
Radius = Radius * 0.001;
moi_tensor = [ m/4*(1/3*h^2+Radius^2) , m/4*(1/3*h^2+Radius^2), m/2*Radius^2];