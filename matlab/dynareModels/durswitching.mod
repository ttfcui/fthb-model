var c d lambda ucost x y a;
varexo eps_y;
varexo_det pshk rshk thetashk deltashk;
parameters beta delta R theta p psi gamma rho mean_y sigma_y;

beta = 0.95;
psi = 1.0;
delta = 0.73;
gamma = 2.0;
R = 1.02;
theta = .20;
p = 1.0;
rho = 0.95;
sigma_y = 0.10;
mean_y = 0.0;

model;

1/exp(c(+1))^(gamma) = (1/(exp(c)^(gamma)) - lambda)/(beta*exp(rshk)*R);
ucost = (exp(rshk)*R-1 + exp(deltashk)*delta)/(exp(rshk)*R);
psi/(exp(pshk)*p*exp(d)^(gamma)) = 
    ucost/(exp(c)^(gamma)) + (exp(thetashk)*theta - ucost)*lambda;
lambda*(exp(x) - exp(c) -
    ((1-exp(deltashk)*delta)/(exp(rshk)*R))*exp(thetashk)*theta*exp(pshk)*p*exp(d)) = 0;
a = exp(x) - exp(c) - exp(pshk)*p*exp(d);
%exp(x(+1)) = exp(rshk)*R*(exp(x) - exp(c)) + 
%    ucost*exp(rshk)*R*exp(pshk)*p*exp(d) + exp(y(+1));
exp(x) = exp(rshk(-1))*R*(exp(x(-1)) - exp(c(-1))) + 
    ucost(-1)*exp(rshk(-1))*R*exp(pshk(-1))*p*exp(d(-1)) + exp(y);
y = (1-rho)*mean_y + rho*y(-1) + eps_y;

end;

initval;

c = 0.0;
d = 0.0;
lambda = -0.1;
ucost = 0.05;
a = 0;
x = 0.0;
y = 0.0;


end;
steady(maxit=5000);
check;

shocks;

var eps_y = sigma_y^2;
var pshk;
periods 1:1;
values -0.00;
var rshk;
periods 1:1;
values -0.10;
var thetashk;
periods 1:1;
values -0.00;
var deltashk;
periods 1:1;
values 0.00;
end;

stoch_simul(irf=0);

time=15;
forecast(periods=15, conf_sig = 0.683) c d lambda a;

oo_.forecast.Mean.c = oo_.forecast.Mean.c - oo_.steady_state(1);
oo_.forecast.HPDsup.c = oo_.forecast.HPDsup.c - oo_.steady_state(1);
oo_.forecast.HPDinf.c = oo_.forecast.HPDinf.c - oo_.steady_state(1);

oo_.forecast.Mean.d = oo_.forecast.Mean.d - oo_.steady_state(2);
oo_.forecast.HPDsup.d = oo_.forecast.HPDsup.d - oo_.steady_state(2);
oo_.forecast.HPDinf.d = oo_.forecast.HPDinf.d - oo_.steady_state(2);

oo_.forecast.Mean.l_fo = ...
    oo_.forecast.Mean.lambda(1:time-1).*exp(oo_.forecast.Mean.c(2:time)).^(gamma);

figure;
subplot(1,2,1);
plot(1:time, [oo_.forecast.Mean.c, ...
     oo_.forecast.HPDsup.c, oo_.forecast.HPDinf.c], ['b']);
title('Forecast For Nondurables');
[oo_.forecast.Mean.c(1:3), oo_.forecast.HPDsup.c(1:3), oo_.forecast.HPDinf.c(1:3)]
oo_.forecast.Mean.l_fo(1:3)
subplot(1,2,2);
plot(1:time, [oo_.forecast.Mean.d, ...
     oo_.forecast.HPDsup.d, oo_.forecast.HPDinf.d], ['b']);
title('Forecast For Durables');
[oo_.forecast.Mean.d(1:3), oo_.forecast.HPDsup.d(1:3), oo_.forecast.HPDinf.d(1:3)]