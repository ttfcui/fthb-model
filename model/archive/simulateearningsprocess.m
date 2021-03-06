% DESCRIPTION:
% This script simulates stochastic income processes, then regresses average
% lifetime earnings on earnings one year before retirement to get (b_0, b_1).
% In the model, predicted lifetime earnings cond. on retirement earnings
% b_0 + b_1*Y_R^i is used to calculate pension schemes as in
% Guvenen and Smith (2014).


%sigmaeta=.01^.5;
%sigmainitial=.15^.5;
%sigmaeps=.05^.5;

%for t=1:35
%    for i=1:100000
%        if t==1
%            z(i,1)=sigmainitial*randn;
%        end
%        z(i,t+1)=z(i,t)+sigmaeta*randn;
%        logearnings(i,t)=z(i,t)+sigmaeps*randn;
%    end
%    sigmaz(t)=std(z(:,t));
%end

%for i=1:100000
%    aveearnings(i)=mean(logearnings(i,:));
%end

%aveearnings=aveearnings';
%retirementearnings=logearnings(:,35);
%y=aveearnings;
%x=[ones(100000,1) retirementearnings];
%beta=inv(x'*x)*(x'*y);
%yhat=x*beta;


% Stochastic component parameters
sigma_z=.21;
rho_z=.91;

% Deterministic component parameters
if exist('Tretire', 'var')
    age = zeros(Tretire, 1);
else
    age = load('../model/input_data/ageearnings.txt');
    Tretire=size(age, 1);
end
n_sim=5e5;
rng(10011)

sigmainitial=((sigma_z^2)/(1-rho_z^2))^.5;
shocks = randn(n_sim,Tretire+1,2);

for t=1:Tretire
    if t==1
        z(:,t) = sigmainitial*shocks(:,t,1);
    end
    z(:,t+1)=rho_z*z(:,t)+sigma_z*shocks(:,t+1,1);
    logearnings(:,t)=z(:,t) + age(t);  % + sigma_eps*shocks(:,t,2);
    sigmaz(t)=std(z(:,t));
end
aveearnings = mean(logearnings, 2); % Age effects get differenced out here
sigmaz

% Note: because shocks to income process are iid & ind. of z_t,
% they do not change the coefficient at all (?)
retirementearnings=logearnings(:,Tretire);
y=aveearnings;
x=[ones(n_sim,1) retirementearnings];
betas=regress(y, x);

% betas(1) can be excluded as it's very close to 0
output = [rho_z, sigma_z, Tretire, betas(2)];
save ../model/input_data/lifeearn.txt output -ASCII;
