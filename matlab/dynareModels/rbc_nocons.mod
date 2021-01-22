// A RBC model with representative agents and firms, and two choice variables:
// Level of durable and level of capital.
// Durables accumulate like capital, but enter in the individual's utility
// function.

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


var k, i, d, x, r, w, y, T, discountK, discountD;
varexo sub; // Percentage subsidy on durable

parameters gamma, beta, alpha, deltaK, delta, taxL, taxK, pD;

gamma  = 3.0000;
beta   = 0.99;
alpha  = 0.36;
deltaK  = 0.05;
delta = deltaK;
taxL   = 0.5;
taxK   = 0.2;
pD = 1.0;

// Third equation in model needs explanation. Find the FOCs for D(+1) and D,
// then use the FOCs on k to reduce the D FOCs into functions of lambda_t.
// Divide the FOCs to get the result.

model(bytecode);
discountK = ((1-taxK)*r+ deltaK*taxK + (1-deltaK));
discountD = pD*(1-sub)*(discountK - (1 - delta));
(d(+1)/d)^(-gamma) = (pD*(1-sub)*(1 - (1 - delta)/discountK))/discountD(-1);
r = alpha*(1/k)^(1-alpha);
w = (1-alpha)*(k/1)^alpha;
x = k - (1-deltaK)*k(-1);
i = d - (1-delta)*d(-1);
y = (k^alpha)*(1^(1-alpha));
T = taxL*w*1 + taxK*r*k - deltaK*taxK*k - sub*i;
x = (1-taxK)*r*k + (1-taxL)*w*1 + deltaK*taxK*k + T - (1-sub)*i;
end;

initval;
sub = 0.01*0.03;
r = (1/beta - deltaK*taxK - (1-deltaK))/(1-taxK);
w = (1-alpha)*(alpha/r)^(alpha/(1-alpha));
discountK = ((1-taxK)*r+ deltaK*taxK + (1-deltaK));
discountD = pD*(1-sub)*(discountK - (1 - delta));
d = 1.5;
T = 0.0;
k = 1.25;
y = 1.0;
x = 0.1;
i = 0.15;
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
periods 9:20;
values 0.3;
end;

simul(periods=250, stack_solve_algo=3);
figure;
hold on
stairs(k(6:105))
stairs(d(6:105), '--', 'LineWidth',0.9)
legend('Capital stock','Durables stock')
title('Model dynamics, without consumption')
print('rbc/Output/RBC_nocons_stock', '-dpdf')
figure;
hold on
stairs(x(6:105))
stairs(i(6:105), '--', 'LineWidth',0.9)
legend('Capital investment','Durables investment')
title('Model dynamics, without consumption')
print('rbc/Output/RBC_nocons_invest', '-dpdf')
