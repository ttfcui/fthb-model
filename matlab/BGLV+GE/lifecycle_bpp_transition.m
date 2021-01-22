function [c_polTrans, h_polTrans, VprimeTrans] = lifecycle_bpp_transition(PVec, VprimeStat)

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

% numerical parameters
tol = 1e-10;

% Transition Years
transYr = length(PVec);
DP = -.15;

% loading income and death probability data
death_data = load('death_data', 'deathprob');
income_data = load('income_data', 'ageearnings');
ageearnings = income_data.ageearnings;

%% Stochastic and deterministic processes

% Functions

function ret_inc = compute_ret_inc
%--------------------------------------------------------------------------
%   Computes retirement income as in the fortran code
%
%   Date created: Dec. 22, 2014
%   Date last modified: Dec. 22, 2014
%   S. Indarte (indarte@u.northwestern.edu)
%--------------------------------------------------------------------------

avg_inc     = 0.0;
pred_inc    = zeros(size(Z));  % pre-allocate
ret_inc     = pred_inc;        % pre-allocate

for j = 1:S
    pred_inc(j)     = 0.3083*(Z(j));
    pred_inc(j)     = exp(pred_inc(j)) / exp(avg_inc);
    
    if (pred_inc(j)<=0.3)
        ret_inc(j)  = 0.9*pred_inc(j)*exp(avg_inc);
    elseif (pred_inc(j)>0.3 && pred_inc(j)<=2) 
        ret_inc(j)  = (0.27+0.32*(pred_inc(j)-0.3))*exp(avg_inc);
    elseif (pred_inc(j)>2 && pred_inc(j)<=4.1) 
        ret_inc(j)  = (0.81+0.15*(pred_inc(j)-2))*exp(avg_inc);
    else
        ret_inc(j)  = 1.13*exp(avg_inc);
    end
end 
    
end



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
    ret_inc = compute_ret_inc;
    y(:,work_years+1:end) = repmat(ret_inc, 1, ret_years+1);
end

% death process
if stoch_death == 1;
    d  = [zeros(1, work_years) death_data.deathprob];
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

Ys = zeros(I,S,lifespan);
Ws = zeros(I,S,lifespan);
for age = 1:lifespan
for s=1:S
    Ys(:,s,age) = y(s,age);
    Ws(:,s,age) = w_grid;
end    
end


%% Policy computation function

VprimeTrans = zeros(I,S,lifespan+1,transYr);
c_polTrans = zeros(I,S,lifespan,transYr);
h_polTrans = c_polTrans;
VprimeTrans(:,:,:,transYr) = VprimeStat;

function [c_polT, h_polT, VprimeT] = policy_transition(P, Vprime, subsidy)

thetaPol = 1 - (1-theta)/exp(DP);
P = exp(DP)*P;    

rent = P*(1-(1-delta)/(1+r));
rent_the = P*(1-(1-thetaPol)*(1-delta)/(1+r));
wprime = w_grid(w_grid>0);
const = (alpha^alpha)*((1-alpha)^(1-alpha));

%--------------------------------------------------------------------------
%   Solve for policy functions at constant prices
%--------------------------------------------------------------------------

Vprime_end = Vprime(:,:,lifespan+1);

% Iterate on age
for age = lifespan:-1:1;
    Vp = Vprime(:,:,age+1);
    if age < work_years
        EVp = Vp*Pr';
    else
        EVp = (1-d(age))*Vp + d(age)*Vprime_end;
    end
    EVp = EVp(w_grid>0,:);
    for s = 1:S;
        % unconstrained
        RHS = const^(sigma-1)*rent^((1-alpha)*(1-sigma))*beta*(1+r)*EVp(:,s);
        c = alpha*(RHS.^(-1/sigma));
        h = (1-alpha)/rent*(RHS.^(-1/sigma));
        check = wprime - thetaPol*(1-delta)*h;
        
        % constrained
        idx = find(check < 0); % indices for which constraint binds            
        h(idx) = wprime(idx)/(thetaPol*(1-delta));
        c_l = alpha/(1-alpha)*rent*h(idx);
        c_h = c(idx);
        F = @(C) ((1-alpha)*h(idx).^((1-alpha)*(1-sigma)-1).*C.^(alpha*(1-sigma))-...
                rent_the*alpha*h(idx).^((1-alpha)*(1-sigma)).*C.^(alpha*(1-sigma)-1))+...
                thetaPol*(1-delta)*beta*EVp(idx,s);
        dif = 1;
        while dif>tol
            c_c = (c_l+c_h)/2;
            F_eval = F(c_c);    
            c_h(F_eval>=0)=c_c(F_eval>=0);
            c_l(F_eval<0)=c_c(F_eval<0);
            dif = max(abs(c_h-c_l));
        end                  
        c(idx) = c_c;

    
        w = wprime/(1+r)+c+rent*h-y(s,age);
        winv(w_grid>0,s,age) = w; 
        winv(w_grid<=0,s,age) = min(w(1),wl);        
        c_polT(:,s,age) = interp1(w,c,w_grid,'linear','extrap');
        h_polT(:,s,age) = interp1(w,h,w_grid,'linear','extrap');
        VprimeT(:,s,age) = alpha*h_polT(:,s,age).^((1-alpha)*(1-sigma)).*...
            c_polT(:,s,age).^(alpha*(1-sigma)-1);
    end
end

end

% lump sum transfer to replicate retirement assets
y(:,work_years+1) = ret_wealth*y(:,work_years)+y(:,work_years+1);  

% Iteration over transition years
subVec = [0, DP, zeros(1, transYr-3)];
for yrs=transYr-1:-1:1
    display(yrs)
    [c_polT, h_polT, V_primeT] = policy_transition(PVec(yrs), ...
        VprimeTrans(:,:,:,yrs+1), subVec(yrs));
     c_polTrans(:,:,:,yrs) = c_polT;
     h_polTrans(:,:,:,yrs) = h_polT;
     VprimeTrans(:,:,1:lifespan,yrs) = V_primeT;
     VprimeTrans(:,:,lifespan+1,yrs) = VprimeStat(:,:,lifespan+1);
end

end