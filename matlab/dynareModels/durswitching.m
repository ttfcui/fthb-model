%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'durswitching';
M_.dynare_version = '4.5.6';
oo_.dynare_version = '4.5.6';
options_.dynare_version = '4.5.6';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('durswitching.log');
M_.exo_names = 'eps_y';
M_.exo_names_tex = 'eps\_y';
M_.exo_names_long = 'eps_y';
M_.exo_det_names = 'pshk';
M_.exo_det_names_tex = 'pshk';
M_.exo_det_names_long = 'pshk';
M_.exo_det_names = char(M_.exo_det_names, 'rshk');
M_.exo_det_names_tex = char(M_.exo_det_names_tex, 'rshk');
M_.exo_det_names_long = char(M_.exo_det_names_long, 'rshk');
M_.exo_det_names = char(M_.exo_det_names, 'thetashk');
M_.exo_det_names_tex = char(M_.exo_det_names_tex, 'thetashk');
M_.exo_det_names_long = char(M_.exo_det_names_long, 'thetashk');
M_.exo_det_names = char(M_.exo_det_names, 'deltashk');
M_.exo_det_names_tex = char(M_.exo_det_names_tex, 'deltashk');
M_.exo_det_names_long = char(M_.exo_det_names_long, 'deltashk');
M_.exo_det_partitions = struct();
M_.endo_names = 'c';
M_.endo_names_tex = 'c';
M_.endo_names_long = 'c';
M_.endo_names = char(M_.endo_names, 'd');
M_.endo_names_tex = char(M_.endo_names_tex, 'd');
M_.endo_names_long = char(M_.endo_names_long, 'd');
M_.endo_names = char(M_.endo_names, 'lambda');
M_.endo_names_tex = char(M_.endo_names_tex, 'lambda');
M_.endo_names_long = char(M_.endo_names_long, 'lambda');
M_.endo_names = char(M_.endo_names, 'ucost');
M_.endo_names_tex = char(M_.endo_names_tex, 'ucost');
M_.endo_names_long = char(M_.endo_names_long, 'ucost');
M_.endo_names = char(M_.endo_names, 'x');
M_.endo_names_tex = char(M_.endo_names_tex, 'x');
M_.endo_names_long = char(M_.endo_names_long, 'x');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names_long = char(M_.endo_names_long, 'a');
M_.endo_partitions = struct();
M_.param_names = 'beta';
M_.param_names_tex = 'beta';
M_.param_names_long = 'beta';
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'R');
M_.param_names_tex = char(M_.param_names_tex, 'R');
M_.param_names_long = char(M_.param_names_long, 'R');
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, 'theta');
M_.param_names_long = char(M_.param_names_long, 'theta');
M_.param_names = char(M_.param_names, 'p');
M_.param_names_tex = char(M_.param_names_tex, 'p');
M_.param_names_long = char(M_.param_names_long, 'p');
M_.param_names = char(M_.param_names, 'psi');
M_.param_names_tex = char(M_.param_names_tex, 'psi');
M_.param_names_long = char(M_.param_names_long, 'psi');
M_.param_names = char(M_.param_names, 'gamma');
M_.param_names_tex = char(M_.param_names_tex, 'gamma');
M_.param_names_long = char(M_.param_names_long, 'gamma');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'mean_y');
M_.param_names_tex = char(M_.param_names_tex, 'mean\_y');
M_.param_names_long = char(M_.param_names_long, 'mean_y');
M_.param_names = char(M_.param_names, 'sigma_y');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_y');
M_.param_names_long = char(M_.param_names_long, 'sigma_y');
M_.param_partitions = struct();
M_.exo_det_nbr = 4;
M_.exo_nbr = 1;
M_.endo_nbr = 7;
M_.param_nbr = 10;
M_.orig_endo_nbr = 7;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 0;
erase_compiled_function('durswitching_static');
erase_compiled_function('durswitching_dynamic');
M_.orig_eq_nbr = 7;
M_.eq_nbr = 7;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 1 6 13;
 2 7 0;
 0 8 0;
 3 9 0;
 4 10 0;
 5 11 0;
 0 12 0;]';
M_.nstatic = 2;
M_.nfwrd   = 0;
M_.npred   = 4;
M_.nboth   = 1;
M_.nsfwrd   = 1;
M_.nspred   = 5;
M_.ndynamic   = 5;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(7, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.maximum_exo_det_lag = 1;
M_.maximum_exo_det_lead = 0;
oo_.exo_det_steady_state = zeros(4, 1);
M_.params = NaN(10, 1);
M_.NNZDerivatives = [37; 93; -1];
M_.params( 1 ) = 0.95;
beta = M_.params( 1 );
M_.params( 6 ) = 1.0;
psi = M_.params( 6 );
M_.params( 2 ) = 0.73;
delta = M_.params( 2 );
M_.params( 7 ) = 2.0;
gamma = M_.params( 7 );
M_.params( 3 ) = 1.02;
R = M_.params( 3 );
M_.params( 4 ) = .20;
theta = M_.params( 4 );
M_.params( 5 ) = 1.0;
p = M_.params( 5 );
M_.params( 8 ) = 0.95;
rho = M_.params( 8 );
M_.params( 10 ) = 0.10;
sigma_y = M_.params( 10 );
M_.params( 9 ) = 0.0;
mean_y = M_.params( 9 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = 0.0;
oo_.steady_state( 2 ) = 0.0;
oo_.steady_state( 3 ) = (-0.1);
oo_.steady_state( 4 ) = 0.05;
oo_.steady_state( 7 ) = 0;
oo_.steady_state( 5 ) = 0.0;
oo_.steady_state( 6 ) = 0.0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
options_.steady.maxit = 5000;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.det_shocks = [ M_.det_shocks;
struct('exo_det',1,'exo_id',1,'multiplicative',0,'periods',1:1,'value',(-0.00)) ];
M_.det_shocks = [ M_.det_shocks;
struct('exo_det',1,'exo_id',2,'multiplicative',0,'periods',1:1,'value',(-0.10)) ];
M_.det_shocks = [ M_.det_shocks;
struct('exo_det',1,'exo_id',3,'multiplicative',0,'periods',1:1,'value',(-0.00)) ];
M_.det_shocks = [ M_.det_shocks;
struct('exo_det',1,'exo_id',4,'multiplicative',0,'periods',1:1,'value',0.00) ];
M_.exo_det_length = 1;
M_.Sigma_e(1, 1) = M_.params(10)^2;
options_.irf = 0;
var_list_ = char();
info = stoch_simul(var_list_);
time=15;
options_.forecasts.conf_sig = 0.683;
options_.periods = 15;
var_list_ = char('c','d','lambda','a');
[oo_.forecast,info] = dyn_forecast(var_list_,M_,options_,oo_,'simul');
oo_.forecast.Mean.c = oo_.forecast.Mean.c - oo_.steady_state(1);
oo_.forecast.HPDsup.c = oo_.forecast.HPDsup.c - oo_.steady_state(1);
oo_.forecast.HPDinf.c = oo_.forecast.HPDinf.c - oo_.steady_state(1);
oo_.forecast.Mean.d = oo_.forecast.Mean.d - oo_.steady_state(2);
oo_.forecast.HPDsup.d = oo_.forecast.HPDsup.d - oo_.steady_state(2);
oo_.forecast.HPDinf.d = oo_.forecast.HPDinf.d - oo_.steady_state(2);
oo_.forecast.Mean.l_fo = ...
    oo_.forecast.Mean.lambda(1:time-1).*exp(oo_.forecast.Mean.c(2:time)).^(gamma);
figure;
subplot(1,2,1);
plot(1:time, [oo_.forecast.Mean.c, ...
     oo_.forecast.HPDsup.c, oo_.forecast.HPDinf.c], ['b']);
title('Forecast For Nondurables');
[oo_.forecast.Mean.c(1:3), oo_.forecast.HPDsup.c(1:3), oo_.forecast.HPDinf.c(1:3)]
oo_.forecast.Mean.l_fo(1:3)
subplot(1,2,2);
plot(1:time, [oo_.forecast.Mean.d, ...
     oo_.forecast.HPDsup.d, oo_.forecast.HPDinf.d], ['b']);
title('Forecast For Durables');
[oo_.forecast.Mean.d(1:3), oo_.forecast.HPDsup.d(1:3), oo_.forecast.HPDinf.d(1:3)]
save('durswitching_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('durswitching_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('durswitching_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('durswitching_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('durswitching_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('durswitching_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('durswitching_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
