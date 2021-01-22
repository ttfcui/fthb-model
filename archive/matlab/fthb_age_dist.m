% This program computes the age distribution of first time home buyers. The
% model is explained in the accompanied PDF.
clear
close all
clc

p_H = 10*exp(1);                        % House price
theta = 0.25;                           % Down payment requirement
mu = 1/19;                              % Drift on log wealth
sigma = 0.0521;                           % Multiplier on BM

a = log(theta*p_H)/sigma;
c = mu/sigma;

t = linspace(10^-4, 45, 1000);          % Time grid
f_tau = a./((2*pi*t.^3).^(1/2)) .* ...  % Density
        exp(-(a-c*t).^2./(2*t));
    
figure;
plot(20+t, f_tau, 'linewidth', 2);
xlabel('Age', 'FontSize', 18)
xlim(20 + [t(1), t(end)])
ylabel('Density', 'FontSize', 18)
set(gca, 'fontsize', 16)
set(gcf, 'paperpositionmode', 'auto')
print('-depsc', 'fthb_age_dist')

[~, max_ind] = max(f_tau);
disp(['The mode is ' num2str(20+t(max_ind))])
disp(['The mean is ' num2str(t*f_tau'/sum(f_tau) + 20)])


%% Comparatice statics

% p_H
disp(' ')
disp('p_H = 15e:')
p_H_cs = 10.1*exp(1);

a_cs = log(theta*p_H_cs)/sigma;

f_tau1 = a_cs./((2*pi*t.^3).^(1/2)) .* ...  % Density
        exp(-(a_cs-c*t).^2./(2*t));
    
figure;
plot(20+t, f_tau, '--', 20+t, f_tau1, 'linewidth', 2);
xlim(20 + [t(1), t(end)])
ylabel('Density', 'FontSize', 18)
set(gca, 'fontsize', 16)
set(gcf, 'paperpositionmode', 'auto')
print('-depsc', 'CS1')

[~, max_ind] = max(f_tau1);
disp(['The mode is ' num2str(20+t(max_ind))])
disp(['The mean is ' num2str(t*f_tau1'/sum(f_tau1) + 20)])

    
    
% theta
disp(' ')
disp('theta')
theta_cs = 0.15;

a_cs = log(theta_cs*p_H)/sigma;

f_tau1 = a_cs./((2*pi*t.^3).^(1/2)) .* ...  % Density
        exp(-(a_cs-c*t).^2./(2*t));
    
figure;
plot(20+t, f_tau, '--', 20+t, f_tau1, 'linewidth', 2);
xlim(20 + [t(1), t(end)])
set(gca, 'fontsize', 16)
set(gcf, 'paperpositionmode', 'auto')
print('-depsc', 'CS2')

[~, max_ind] = max(f_tau1);
disp(['The mode is ' num2str(20+t(max_ind))])
disp(['The mean is ' num2str(t*f_tau1'/sum(f_tau1) + 20)])
    
% mu
disp(' ')
disp('mu = 1/15')
mu_cs = 1/15;

c_cs = mu_cs/sigma;

f_tau1 = a./((2*pi*t.^3).^(1/2)) .* ...  % Density
        exp(-(a-c_cs*t).^2./(2*t));
    
figure;
plot(20+t, f_tau, '--', 20+t, f_tau1, 'linewidth', 2);
xlabel('Age', 'FontSize', 18)
xlim(20 + [t(1), t(end)])
ylabel('Density', 'FontSize', 18)
set(gca, 'fontsize', 16)
set(gcf, 'paperpositionmode', 'auto')
print('-depsc', 'CS3')

[~, max_ind] = max(f_tau1);
disp(['The mode is ' num2str(20+t(max_ind))])
disp(['The mean is ' num2str(t*f_tau1'/sum(f_tau1) + 20)])

% sigma
disp(' ')
disp('sigma = 0.2')
sigma_cs = 0.2;

a_cs = log(theta*p_H)/sigma_cs;
c_cs = mu/sigma_cs;

f_tau1 = a_cs./((2*pi*t.^3).^(1/2)) .* ...  % Density
        exp(-(a_cs-c_cs*t).^2./(2*t));
    
figure;
plot(20+t, f_tau, '--', 20+t, f_tau1, 'linewidth', 2);
xlabel('Age', 'FontSize', 18)
xlim(20 + [t(1), t(end)])
set(gca, 'fontsize', 16)
set(gcf, 'paperpositionmode', 'auto')
print('-depsc', 'CS4')

[~, max_ind] = max(f_tau1);
disp(['The mode is ' num2str(20+t(max_ind))])
disp(['The mean is ' num2str(t*f_tau1'/sum(f_tau1) + 20)])

%% Brownian Motion Example

dt = 1/1000;
T = 40;
w = zeros(1+T/dt, 1);
dB_t = sigma*sqrt(dt)*randn(T/dt, 1);

for i = 1:T/dt
    w(i+1) = w(i) + mu*dt + dB_t(i);
end

figure; hold on
plot(20+ [0:dt:T], w)
plot(20+ [0:dt:T], ones(1, T/dt+1)*log(p_H*theta), '--k')