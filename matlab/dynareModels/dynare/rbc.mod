// A RBC model with representative agents and firms, and three choice variables:
// Level of consumption, level of capital and level of durable.
// Durables enter the utility function only. Capital enters the production function only.
// A preliminary housing market clearing condition is included (last equation).


/*
 * Copyright (C) 2001-2010 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */


var c, k, i, d, x, r, w, y, T, pD, discountK;
varexo sub; // Percentage subsidy on durable.

parameters gamma, beta, alpha, alphaU, deltaK, delta, taxL, taxK, psi1, psi2;

gamma  = 3.0;
beta   = 0.98;
alpha  = 0.36;
alphaU = 0.80;
deltaK  = 0.05;
delta = 0.05;
taxL   = 0.3;
taxK   = 0.2;
psi1 = 1.0;
psi2 = 2.5/(1+2.5);


model(bytecode);
discountK = ((1-taxK)*r(+1)+ (1-deltaK));
d(+1) = ((1-alphaU)/alphaU)*c(+1)/
    (pD*(1-sub)*discountK - pD(+1)*(1-sub(+1))*(1-delta));
(d/c) = (beta*((c(+1)^(alphaU)*d(+1)^(1-alphaU))/
    (c^(alphaU)*d^(1-alphaU)))^(-gamma)*
    discountK)^(1/(1-alphaU))*(d(+1)/c(+1));
r = alpha*(1/k)^(1-alpha);
w = (1-alpha)*(k/1)^alpha;
x = k - (1-deltaK)*k(-1);
i = d - (1-delta)*d(-1);
y = (k^alpha)*(1^(1-alpha));
T = taxL*w*1 + taxK*r*k - pD*sub*i;
c = (1-taxK)*r*k + (1-taxL)*w*1 + T - x - pD*(1-sub)*i;
d = (1-delta)*d(-1) + psi1*(pD*psi1*psi2)^(-psi2/(psi2-1));
end;

initval;
sub = 0.01*0.03;
r = (1/beta - deltaK*taxK - (1-deltaK))/(1-taxK);
w = (1-alpha)*(alpha/r)^(alpha/(1-alpha));
discountK = ((1-taxK)*r+ deltaK*taxK + (1-deltaK));
c = 0.2;
d = 1.5;
T = 0.0;
k = 1.25;
y = 1.0;
x = 0.1;
i = 0.15;
pD = 1.1;
end;
steady(maxit=500);
check;

// endval;
// sub = 0.015;
// c = 0.2;
// k = 0.5;
// r = ((1-sub)/beta - (delta*taxK + (1-sub)*(1-delta)))/(1-taxK);
// w = (1-alpha)*(alpha/r)^(alpha/(1-alpha));
// x = 0.3;
// end;
// steady;

shocks;
var sub;
periods 90:90;
values 0.06;
end;

simul(periods=200, stack_solve_algo=3);
// Export graphical output
figure;
hold on
stairs(d(6:150), '--', 'LineWidth',0.9)
legend('Durables stock')
title('Model dynamics, with consumption')
print('rbc/Output/RBC_stock', '-dpdf')
figure;
hold on
stairs(i(6:150),'--', 'LineWidth',0.9)
legend('Durables investment')
title('Model dynamics, with consumption')
print('rbc/Output/RBC_invest', '-dpdf')
figure;
hold on
stairs(c(6:150), '-k')
legend('Consumption dynamics')
print('rbc/Output/RBC_cons', '-dpdf')
