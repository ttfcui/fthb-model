% Reproduces age distributions from DB's old script
clear
close all
clc

orig = fthb_age_dist(10*exp(1), 0.1, 1/12, 0.15, 'fthb_age_dist');

%% Comparative statics

% p_H
disp(' ')
disp('p_H = 15e:')
fthb_age_dist(15*exp(1), 0.1, 1/12, 0.15, 'CS1', orig);
        
% theta
disp(' ')
disp('theta = 0.15:')
fthb_age_dist(10*exp(1), 0.15, 1/12, 0.15, 'CS2', orig);

% mu
disp(' ')
disp('mu = 1/15')
fthb_age_dist(10*exp(1), 0.1, 1/15, 0.15, 'CS3', orig);

% sigma
disp(' ')
disp('sigma = 0.2')
fthb_age_dist(10*exp(1), 0.1, 1/12, 0.2, 'CS4', orig);

%% Wider range of CS for just p_H

for i = 5:15
    disp(' ')
    fthb_age_dist(i*exp(1), 0.1, 1/12, 0.15, strcat('CS_pH_',num2str(i)), ...
                  orig);
    disp(' ')
    fthb_age_dist(10*exp(1), 0.1, 1/i, 0.15, strcat('CS_mu_',num2str(i)), ...
                  orig);
end


%% One simulation of a Brownian Motion
%% (The dotted line is the cost of a house)

fthb_age_simul(10*exp(1), 0.1, 1/12, 0.15)
