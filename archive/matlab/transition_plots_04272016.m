% This script loads the simulation results from fortran and plots the
% fthb steady state distribution and the transition dynamics
close all
clc
clear 


fthb_ss = importdata('dist_fthb.txt');
hstock_orig = importdata('housingstock.txt');
T = 40;                         % number of periods
fthb_ss = fthb_ss(fthb_ss(:,2) <= T,:);
hstock_orig = hstock_orig(1:T,:);



%% Steady State Plots
H = size(fthb_ss, 1);           % number of HH who own a house during their lifetime
h_ss = sum(hstock_orig(:, 2)); % Total housing wealth in steady state
i_ss = sum(hstock_orig(:, 3)); % Total housing *investment* in steady state
q_ss = H;               % number of first time houses bought in every period
age = 19+ [1:T];
fthb_age = zeros(T, 1);
for t = 1:T
    fthb_age(t) = sum(fthb_ss(:, 2)==t)/H;
end

% Age distribution
figure;
plot(age, fthb_age, 'linewidth', 2)
set(gca, 'fontsize', 16)
title('Age Distribution')
ylabel('% of FTH bought by age')
xlabel('Age')
xlim([age(1), age(end)])
print('fracFthb','-dpdf')

% Size distribution
figure;
hist(fthb_ss(:, 3), 50)
xlabel('Size')
title('Histogram of fth size')
set(gca, 'fontsize', 16)
print('HouseSize','-dpdf')

% Wealth distribution
figure;
hist(fthb_ss(:, 4), 50)
xlabel('Size')
title('Histogram of fthb wealth')
set(gca, 'fontsize', 16)
print('WealthHist', '-dpdf')

% Income distribution
figure;
hist(fthb_ss(:, 5), 50)
xlabel('Size')
title('Histogram of fthb income')
set(gca, 'fontsize', 16)
print('IncHist', '-dpdf')

% 3D Scatter of age, wealth + income, and size when fth bought
% figure;
% scatter3(19 + fthb_ss(:, 2), fthb_ss(:, 3), fthb_ss(:, 4) + fthb_ss(:, 5))
% xlabel('Age')
% ylabel('Size')
% zlabel('Wealth + income')
% title('Scatter of FTHB')
% set(gca, 'fontsize', 16)
% print('AssetScatter', '-dpdf')

%% Transition to a policy change

fthb_transition =  importdata('transition_fthb.txt');
hstock = importdata('housing_transit.txt');

T_plot = T + 2;
for t = 0:-1: -T_plot
    fthb_append = [fthb_ss(:, 1), ones(H, 1)*t, fthb_ss(:, 2:end)];
    fthb_transition = [fthb_transition; fthb_append];
end

%% Affix all the period-cohort statistics together
j = 4;
PolYrs = 1;
Q_p = zeros(T_plot+1, 1); % FTHB number holder
H_p = zeros(T, 1); % Housing wealth holder
I_p = zeros(T, 1); % Housing investment holder
fthb_age_transition = zeros(T, 2);
for p = 0:T_plot
    
    p_ind = find((fthb_transition(:, 3) - fthb_transition(:, 2)) == p);
    fthb_p = fthb_transition(p_ind, :);
    Q_p(p+1) = size(fthb_p,1)/q_ss;

    if (p < T)
    if (p > 0)
    h_ind = find((hstock(:,1) == p) & (hstock(:,2) > T - p) & (hstock(:,2) <= T));
    hstock(h_ind, 3) = hstock_orig(1:p, 2);     
    hstock(h_ind, 4) = hstock_orig(1:p, 3);     
    end
    H_p(p+1) = sum(hstock((hstock(:, 1) == p), 3))/h_ss;
    I_p(p+1) = sum(hstock((hstock(:, 1) == p), 4))/i_ss;
    end
    H_pd = diff([1; H_p]) + 1; % Delta Housing wealth
    I_pd = diff([1; I_p]) + 1; % Delta Housing investment


end

% Get age distribution over time
for t = 1:T
    p_ind = find(((fthb_transition(:, 3) - fthb_transition(:, 2)) < PolYrs) ...
        & (fthb_transition(:, 3) == t));
    fthb_age_transition(t, 1) = size(fthb_transition(p_ind, 1), 1) ...
        /sum(Q_p(1:PolYrs)*q_ss);
    p_ind = find(((fthb_transition(:, 3) - fthb_transition(:, 2)) >= PolYrs) ...
        & (fthb_transition(:, 3) == t));
    fthb_age_transition(t, 2) = size(fthb_transition(p_ind, 1), 1); % ...
    %    /sum(Q_p(PolYrs+1:T_plot)*q_ss);
end

fthb_age_transition(:, 2) = ...
    fthb_age_transition(:, 2)/sum(fthb_age_transition(:, 2));

figure;
hold on
title('Transaction following a policy change')
stairs(-3:T_plot, [ones(3, 1); Q_p])
stairs(-3:T-1, [ones(3, 1); I_pd], 'r')
plot(-3:T_plot, 1, '-.k')
ylabel('Index normalized, 1 = steady state value')
xlabel('Period')
legend('FTH purchases', 'Growth, housing investment')

lims = get(gca, 'YLim');
if (lims(1) > 0.9)
    lims(1) = 0.9;
end
if (lims(2) < 1.1)
    lims(2) = 1.1;
end

axis([get(gca, 'XLim'), lims])
print('FthbShock','-dpdf')


% Cumulative gains over steady-state flow
Q_pc = cumsum(Q_p - 1);
I_pc = I_p - 1;
figure;
hold on
title('Transaction following a policy change')
stairs(-3:T_plot, [zeros(3, 1); Q_pc])
stairs(-3:T-1, [zeros(3, 1); I_pc], 'r')
plot(-3:T_plot, 0, '-.k')
ylabel('Index normalized, 1 = steady state value')
xlabel('Period')
legend('FTH purchases', 'Growth, housing investment')

lims = get(gca, 'YLim');
if (lims(1) > -0.1)
    lims(1) = -0.1;
end
if (lims(2) < 0.1)
    lims(2) = 0.1;
end

axis([get(gca, 'XLim'), lims])
print('FthbShockCum','-dpdf')


% Stacked bar plot of age distribution
figure;
    subplot(2,1,1);
    hold on
    bar(age+0.5, fthb_age_transition(:,1))
    stairs(age, fthb_age, '-.', 'Color', 'r', 'LineWidth', 1.9)
    eval(['title(''Policy period( to period ' num2str(PolYrs) ')'');'])
    ylim([0 0.1])
    subplot(2,1,2);
    hold on
    bar(age+0.5, fthb_age_transition(:,2))
    stairs(age, fthb_age, '-.', 'Color', 'r', 'LineWidth', 1.9)
    title('Remaining Periods')
    xlabel('Age')
    ylabel('Prop. of FTH bought by age')
    ylim([0 0.1])
%     set(gca, 'fontsize', 16)
print('FthbShockAge','-dpdf')
% plot(age, fthb_age, '--', age, fthb_age_transition(:, 1), age, fthb_age_transition(:, 5), 'linewidth', 3)
% legend('Steady State', 'On impact', 'After 4 periods')
