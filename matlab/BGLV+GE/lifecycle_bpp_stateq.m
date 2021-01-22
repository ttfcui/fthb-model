function [diffH, H_sims, Vprime, probs] = lifecycle_bpp_stateq(P)

%-------------------------------------------------------------------------
%   Life cycle model
%   Compute elasticities, MPC*H and decomposition
%
%   Calibration to match life cycle profile of housing and non-housing
%   wealth (SCF 2001)
%
%   sigma = 2
%
%   Date created: 5/27/2015
%   Note: requires policy_constant_P.m
%-------------------------------------------------------------------------

close all


%% Parameters
% Computational parameters
stoch_death     = 0;    % = 1 if death is to be stochastic
comp_ret_inc    = 1;    % = 1 to compute retirement income as in fortran code

% Model parameters
r           = 0.024;
sigma       = 2; 
delta       = 0.022;
theta       = 0.25;
alpha       = 0.8539;
beta        = 0.9331;
psi         = 1.3529e+3;
ret_wealth  = 3.4137;
work_years  = 35;
ret_years   = 25;
lifespan    = work_years + ret_years;
bins        = 1:5:51;
n_bins      = length(bins);


elas = 2.5;   % Price elasticity of construction
psi = elas/(1+elas);  % elas = psi/(1-psi) in Stroebel paper 

% numerical parameters
tol = 1e-10;
tolP = 1e-5;

% loading income and death probability data
load death_data             % loads a 1 x 25 vector of death probabilities
load income_data            % loads a 1 x 35 vector of age-specific income components
load targets_by_age         % loads a vector of moments from SCF 2001 (home owners)

%% Stochastic and deterministic processes
% Income process
chi                     = zeros(1,lifespan+1); % age-specific component of income
chi(1:work_years)       = ageearnings;
chi(work_years+1:end)   = 0; 

rho_z = 0.91;
sigma_z = 0.21;

[Z, Pr]             = tauchen(13,0,rho_z,sigma_z,2.5);
S                   = length(Z);
pr                  = ones(1,S)/S; dis = 1;
while dis>tol
    pri=pr*Pr; dis=max(abs(pri-pr)); pr=pri;
end
TM = ((pr'*ones(1,S))./(ones(S,1)*pr).*Pr)';

y = zeros(S,lifespan+1);
for i = 1:S
    y(i,:) = exp(Z(i) + chi)/(pr*exp(Z));
end

if comp_ret_inc == 1
    compute_ret_inc
    y(:,work_years+1:end) = repmat(ret_inc, 1, ret_years+1);
end

% death process
if stoch_death == 1;
    d  = [zeros(1, work_years) deathprob];
else
    d = zeros(1,lifespan);
end

% mean and median earnings
Y = pr*y;

ys = y(:,1:work_years);
pr_income = pr'*ones(1,work_years)/work_years;
dist_income = [ys(:) pr_income(:)];
dist_income = sortrows(dist_income);
dist_income(:,2) = cumsum(dist_income(:,2));
med_ear = dist_income(find(dist_income(:,2)>0.5,1),1);

% grids
I = 800; wl = -.1; wh = 30;
w_grid = linspace(wl,wh,I)';
initial_wealth

Ys = zeros(I,S,lifespan);
Ws = zeros(I,S,lifespan);
for age = 1:lifespan
for s=1:S
    Ys(:,s,age) = y(s,age);
    Ws(:,s,age) = w_grid;
end    
end


%% Compute policy functions

% lump sum transfer to replicate retirement assets
y(:,work_years+1) = ret_wealth*y(:,work_years)+y(:,work_years+1);

rent = P*(1-(1-delta)/(1+r));
rent_the = P*(1-(1-theta)*(1-delta)/(1+r));
wprime = w_grid(w_grid>0);
const = (alpha^alpha)*((1-alpha)^(1-alpha));

policy_constant_P



%% Wealth distributions by income and age

M = zeros(I,S,lifespan);
M(w_grid>=0,:,:)=1;
M(:,:,1) = M0;
Mi = zeros(I, S); % pre-allocate

for age = 1:lifespan-1
    for j=1:S
        Mi(:,j) = interp1(w_grid,M(:,j,age),winv(:,j,age),'linear','extrap');
        Mi(:,j) = max(min(Mi(:,j),1),0);
    end
    M(:,:,age+1) = Mi*TM';
end

m = zeros(I,S,lifespan);
m(2:(end),:,:) = M(2:I,:,:) - M(1:I-1,:,:);

% joint distribution of wealth, income shock, age

probs = pr'*ones(1,lifespan)/(lifespan);
probs = permute(repmat(probs,[1 1 I]),[3 1 2]);
probs = m.*probs;

%% Means by age for housing, liquid wealth and retirement wealth

As = Ws + Ys - c_pol - h_pol;
bins        = 1:5:51;
for i_bin = 1:n_bins
    H_bin = h_pol(:,:,bins(i_bin):bins(i_bin)+4);
    A_bin = As(:,:,bins(i_bin):bins(i_bin)+4);
    ps = probs(:,:,bins(i_bin):bins(i_bin)+4);
    EH(i_bin) = H_bin(:)'*ps(:)/sum(ps(:));
    EA(i_bin) = A_bin(:)'*ps(:)/sum(ps(:));
end

load targets_by_age     


bins(9)=mean(bins(9:end));
bins=bins(1:9);
ages = 24+bins;

EH(9)=mean(EH(9:end));
EH=EH(1:9);

EA(9)=mean(EA(9:end));
EA=EA(1:9);






%% Simulations
% Simulation parameters
get_elas = 1;
nPeople = 200;
    % SET INITIAL INCOME
    %H_sims(1,:) = 8;  % was trying to understand what effects MPC
    %W_sims(1,:) = 8;  % it is theta interacted with initial wealth (not
    %housing wealth)
    % cap wealth at 10
C_sims = zeros(work_years,nPeople); 
W_sims = zeros(lifespan+1,nPeople);
H_sims = zeros(work_years,nPeople);
Y_sims = zeros(work_years,nPeople);
inc_idx = zeros(1,nPeople);    

% initialize income state
Probinitcum = cumsum(pr);
sindex = ceil(nPeople*Probinitcum);
Y_sims(1,1:sindex(1)) = y(1,1);
inc_idx(1,1:sindex(1)) = 1;
for s = 2:S
    Y_sims(1,(sindex(s-1)+1:sindex(s))) = y(s,1);
    inc_idx(1,(sindex(s-1)+1:sindex(s))) = s;
end
  
    
   % do simulation 
    for i = 1:nPeople
        for j = 1:lifespan
            k =  inc_idx(i);
            C_sims(j,i) = interp1(w_grid,c_pol(:,k,j), W_sims(j,i),...
                'linear','extrap');
            H_sims(j,i) = interp1(w_grid,h_pol(:,k,j), W_sims(j,i),...
                'linear','extrap');
            if j < work_years % <, not <= b/c it's setting j+1
                [~,kk] = max(rand(1) <= cumsum(Pr(k,:)));
                inc_idx(i) = kk;
                Y_sims(j+1,i) = y(kk, j+1);
            else
                Y_sims(j+1,i) = y(kk, j+1);
            end
        end
        W_sims(j+1,i) = (1+r)*(W_sims(j,i) + Y_sims(j,i) - C_sims(j,i) -...
            rent*H_sims(j,i));
    end
    Y_sims = Y_sims(1:lifespan,:);
    

%% Calculate new housing construction - does it make up for depreciation?

H = sum(sum(H_sims));
diffH = (abs(delta*H - 10*(psi*10*P)^((-psi)/(psi-1))))/H; % psi_1 in Stroebel = 10

end


%%%%%  Now try it for y_sims and c_sims - I think put in companion form and
%%%%%  then run the regression then do it for all the income groups
%%%%%  together and then increase the size of the simulation. I think we
%%%%%  want this to be one big matrix
% set number of firms
