%--------------------------------------------------------------------------
%   Solve for policy functions at constant prices
%--------------------------------------------------------------------------

% Terminal conditions
Vprime_end = psi*const^(1-sigma)*rent^((1-alpha)*(sigma-1))*...
    (ones(I,1)*y(:,end)'/r+w_grid*ones(1,S)).^(-sigma);
Vprime = zeros(I,S,lifespan+1);
Vprime(:,:,lifespan+1) = Vprime_end;

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
        check = wprime - theta*(1-delta)*h;
        
        % constrained
        idx = find(check < 0); % indices for which constraint binds            
        h(idx) = wprime(idx)/(theta*(1-delta));
        c_l = alpha/(1-alpha)*rent*h(idx);
        c_h = c(idx);
        F = @(C) ((1-alpha)*h(idx).^((1-alpha)*(1-sigma)-1).*C.^(alpha*(1-sigma))-...
                rent_the*alpha*h(idx).^((1-alpha)*(1-sigma)).*C.^(alpha*(1-sigma)-1))+...
                theta*(1-delta)*beta*EVp(idx,s);
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
        c_pol(:,s,age) = interp1(w,c,w_grid,'linear','extrap');
        h_pol(:,s,age) = interp1(w,h,w_grid,'linear','extrap');
        Vprime(:,s,age) = alpha*h_pol(:,s,age).^((1-alpha)*(1-sigma)).*...
            c_pol(:,s,age).^(alpha*(1-sigma)-1);
    end
end






