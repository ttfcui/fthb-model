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