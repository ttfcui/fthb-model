% construct distribution of entering young cohort

factor = exp(ageearnings(1));

%only homeowners: (22-24)
%w1 = factor*(0.4003 + theta*2.9545);
%w2 = factor*(1.3628 + theta*1.2727);
%w3 = factor*(2.7340 + theta*4.0909);
%w4 = factor*(0.0909 + theta*3.0455);

%everyone: (22-24)
%w1 = factor*(0.4003 + theta*2.9545)*.0957;
%w2 = factor*((1.3628 + theta*1.2727)*.0883+(1-.0883)*.0023);
%w3 = factor*((2.7340 + theta*4.0909)*.1736+(1-.1736)*.0139);
%w4 = factor*((0.0909 + theta*3.0455)*.2793+(1-.2793)*.0998);

%everyone (23-27):
w1 = factor*(0.3145 + theta*2.3214)*.1156+(1-.1156)*0;
w2 = factor*((0.2556 + theta*1.7857)*.2433+(1-.2433)*.0075);
w3 = factor*((.1611 + theta*2.8571)*.2658+(1-.2658)*.0289);
w4 = factor*((0.4902 + theta*4.0357)*.5406+(1-.5406)*.1689);


i1 = find(w_grid>w1,1);
i2 = find(w_grid>w2,1);
i3 = find(w_grid>w3,1);
i4 = find(w_grid>w4,1);

M0 = zeros(I,S);
M0(i1:end,1:5) = 1;
M0(i2:end,6) = 1;
M0(i2:i3-1,7) = 1/2;
M0(i3:end,7) = 1;
M0(i3:end,8) = 1;
M0(i4:end,9:13) = 1;
