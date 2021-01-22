% This script loads the simulation results from fortran and plots the
% fthb steady state distribution and the transition dynamics
close all
clc
clear


hstock_orig = importdata('housingstock.txt');
T = 40;                         % number of periods
H = sum(hstock_orig(1:T,4));           % number of durables owners in each period
h_ss = sum(hstock_orig(1:T, 2)); % Total housing wealth in steady state
q_ss = H;

hstock = importdata('housing_transit.txt');

%% Affix all the period-cohort statistics together
Q_p = zeros(T, 1); % Total durables holder
H_p = zeros(T, 1); % Housing wealth holder
for p = 0:T

    if (p < T)
    if (p > 0)
    h_ind = find((hstock(:,1) == p) & (hstock(:,2) > T - p) & (hstock(:,2) <= T));
    hstock(h_ind, 3) = hstock_orig(1:p, 2);
    hstock(h_ind, 5) = hstock_orig(1:p, 4);
    end
    H_p(p+1) = sum(hstock((hstock(:, 1) == p), 3))/h_ss;
    Q_p(p+1) = sum(hstock((hstock(:, 1) == p), 5))/q_ss;
    end
    H_pd = diff([1; H_p]) + 1; % Delta Housing wealth

end

% Year-by-year gains
figure;
hold on
title('Transaction following a policy change')
stairs(-3:T-1, [ones(3, 1); Q_p])
stairs(-3:T-1, [ones(3, 1); H_pd], 'r')
plot(-3:T-1, 1, '-.k')
ylabel('Index normalized, 1 = steady state value')
xlabel('Period')
legend('FTH purchases', 'Housing wealth')

lims = get(gca, 'YLim');
if (lims(1) > 0.9)
    lims(1) = 0.9;
end
if (lims(2) < 1.1)
    lims(2) = 1.1;
end

axis([get(gca, 'XLim'), lims])
print('CARSShock','-dpdf')


% Cumulative gains over steady-state flow
Q_pc = cumsum(Q_p - 1);
H_pc = H_p - 1;
figure;
hold on
title('Transaction following a policy change')
stairs(-3:T-1, [zeros(3, 1); Q_pc])
stairs(-3:T-1, [zeros(3, 1); H_pc], 'r')
plot(-3:T-1, 0, '-.k')
ylabel('Index normalized, 1 = steady state value')
xlabel('Period')
legend('FTH purchases', 'Housing wealth')

lims = get(gca, 'YLim');
if (lims(1) > -0.1)
    lims(1) = -0.1;
end
if (lims(2) < 0.1)
    lims(2) = 0.1;
end

axis([get(gca, 'XLim'), lims])
print('CARSShockCum','-dpdf')
