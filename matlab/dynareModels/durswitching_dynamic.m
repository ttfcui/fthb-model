function [residual, g1, g2, g3] = durswitching_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(7, 1);
T11 = exp(y(13))^params(7);
T15 = exp(y(6))^params(7);
T45 = exp(x(it_, 2))*params(5)*exp(y(7))^params(7);
T65 = exp(y(7))*params(5)*exp(x(it_, 2))*params(4)*exp(x(it_, 4))*(1-exp(x(it_, 5))*params(2))/(exp(x(it_, 3))*params(3));
lhs =1/T11;
rhs =(1/T15-y(8))/(params(1)*exp(x(it_, 3))*params(3));
residual(1)= lhs-rhs;
lhs =y(9);
rhs =(exp(x(it_, 3))*params(3)-1+exp(x(it_, 5))*params(2))/(exp(x(it_, 3))*params(3));
residual(2)= lhs-rhs;
lhs =params(6)/T45;
rhs =y(9)/T15+y(8)*(exp(x(it_, 4))*params(4)-y(9));
residual(3)= lhs-rhs;
residual(4) = y(8)*(exp(y(10))-exp(y(6))-T65);
lhs =y(12);
rhs =exp(y(10))-exp(y(6))-exp(x(it_, 2))*params(5)*exp(y(7));
residual(5)= lhs-rhs;
lhs =exp(y(10));
rhs =params(3)*exp(x(it_-1, 3))*(exp(y(4))-exp(y(1)))+params(5)*params(3)*exp(x(it_-1, 3))*y(3)*exp(x(it_-1, 2))*exp(y(2))+exp(y(11));
residual(6)= lhs-rhs;
lhs =y(11);
rhs =(1-params(8))*params(9)+params(8)*y(5)+x(it_, 1);
residual(7)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(7, 18);

  %
  % Jacobian matrix
  %

T111 = exp(y(6))*getPowerDeriv(exp(y(6)),params(7),1);
T124 = exp(y(13))*getPowerDeriv(exp(y(13)),params(7),1);
T128 = (-(params(5)*params(3)*exp(x(it_-1, 3))*y(3)*exp(x(it_-1, 2))*exp(y(2))));
T130 = exp(y(7))*getPowerDeriv(exp(y(7)),params(7),1);
T145 = (-(exp(y(2))*params(5)*params(3)*exp(x(it_-1, 3))*exp(x(it_-1, 2))));
T157 = (-(params(3)*exp(x(it_-1, 3))*(exp(y(4))-exp(y(1)))+params(5)*params(3)*exp(x(it_-1, 3))*y(3)*exp(x(it_-1, 2))*exp(y(2))));
T159 = (-((1/T15-y(8))*params(1)*exp(x(it_, 3))*params(3)));
T160 = params(1)*exp(x(it_, 3))*params(3)*params(1)*exp(x(it_, 3))*params(3);
T163 = exp(x(it_, 3))*params(3)*exp(x(it_, 3))*params(3);
T176 = (-(exp(y(7))*params(5)*exp(x(it_, 2))*params(4)*exp(x(it_, 4))*(-(exp(x(it_, 3))*params(3)*(1-exp(x(it_, 5))*params(2))))/T163));
T177 = y(8)*T176;
T189 = (-(exp(y(7))*params(5)*exp(x(it_, 2))*params(4)*exp(x(it_, 4))*(-(exp(x(it_, 5))*params(2)))/(exp(x(it_, 3))*params(3))));
T190 = y(8)*T189;
  g1(1,6)=(-((-T111)/(T15*T15)/(params(1)*exp(x(it_, 3))*params(3))));
  g1(1,13)=(-T124)/(T11*T11);
  g1(1,8)=(-((-1)/(params(1)*exp(x(it_, 3))*params(3))));
  g1(1,16)=(-(T159/T160));
  g1(2,9)=1;
  g1(2,16)=(-((T163-exp(x(it_, 3))*params(3)*(exp(x(it_, 3))*params(3)-1+exp(x(it_, 5))*params(2)))/T163));
  g1(2,18)=(-(exp(x(it_, 5))*params(2)/(exp(x(it_, 3))*params(3))));
  g1(3,6)=(-((-(y(9)*T111))/(T15*T15)));
  g1(3,7)=(-(params(6)*exp(x(it_, 2))*params(5)*T130))/(T45*T45);
  g1(3,8)=(-(exp(x(it_, 4))*params(4)-y(9)));
  g1(3,9)=(-(1/T15-y(8)));
  g1(3,15)=(-(params(6)*T45))/(T45*T45);
  g1(3,17)=(-(y(8)*exp(x(it_, 4))*params(4)));
  g1(4,6)=y(8)*(-exp(y(6)));
  g1(4,7)=y(8)*(-T65);
  g1(4,8)=exp(y(10))-exp(y(6))-T65;
  g1(4,10)=y(8)*exp(y(10));
  g1(4,15)=y(8)*(-T65);
  g1(4,16)=T177;
  g1(4,17)=y(8)*(-T65);
  g1(4,18)=T190;
  g1(5,6)=exp(y(6));
  g1(5,7)=exp(x(it_, 2))*params(5)*exp(y(7));
  g1(5,10)=(-exp(y(10)));
  g1(5,12)=1;
  g1(5,15)=exp(x(it_, 2))*params(5)*exp(y(7));
  g1(6,1)=(-(params(3)*exp(x(it_-1, 3))*(-exp(y(1)))));
  g1(6,2)=T128;
  g1(6,3)=T145;
  g1(6,4)=(-(params(3)*exp(x(it_-1, 3))*exp(y(4))));
  g1(6,10)=exp(y(10));
  g1(6,11)=(-exp(y(11)));
  g1(6,15)=T128;
  g1(6,16)=T157;
  g1(7,5)=(-params(8));
  g1(7,11)=1;
  g1(7,14)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(93,3);
T194 = T111+exp(y(6))*exp(y(6))*getPowerDeriv(exp(y(6)),params(7),2);
  v2(1,1)=1;
  v2(1,2)=96;
  v2(1,3)=(-((T15*T15*(-T194)-(-T111)*(T15*T111+T15*T111))/(T15*T15*T15*T15)/(params(1)*exp(x(it_, 3))*params(3))));
  v2(2,1)=1;
  v2(2,2)=229;
  v2(2,3)=(T11*T11*(-(T124+exp(y(13))*exp(y(13))*getPowerDeriv(exp(y(13)),params(7),2)))-(-T124)*(T11*T124+T11*T124))/(T11*T11*T11*T11);
  v2(3,1)=1;
  v2(3,2)=276;
  v2(3,3)=(-((-(params(1)*exp(x(it_, 3))*params(3)*(-T111)/(T15*T15)))/T160));
  v2(4,1)=1;
  v2(4,2)=106;
  v2(4,3)=  v2(3,3);
  v2(5,1)=1;
  v2(5,2)=278;
  v2(5,3)=(-(params(1)*exp(x(it_, 3))*params(3)/T160));
  v2(6,1)=1;
  v2(6,2)=142;
  v2(6,3)=  v2(5,3);
  v2(7,1)=1;
  v2(7,2)=286;
  v2(7,3)=(-((T159*T160-T159*(T160+T160))/(T160*T160)));
  v2(8,1)=2;
  v2(8,2)=286;
  v2(8,3)=(-((T163*(T163+T163-(T163+exp(x(it_, 3))*params(3)*(exp(x(it_, 3))*params(3)-1+exp(x(it_, 5))*params(2))))-(T163-exp(x(it_, 3))*params(3)*(exp(x(it_, 3))*params(3)-1+exp(x(it_, 5))*params(2)))*(T163+T163))/(T163*T163)));
  v2(9,1)=2;
  v2(9,2)=322;
  v2(9,3)=(-((-(exp(x(it_, 3))*params(3)*exp(x(it_, 5))*params(2)))/T163));
  v2(10,1)=2;
  v2(10,2)=288;
  v2(10,3)=  v2(9,3);
  v2(11,1)=2;
  v2(11,2)=324;
  v2(11,3)=(-(exp(x(it_, 5))*params(2)/(exp(x(it_, 3))*params(3))));
  v2(12,1)=3;
  v2(12,2)=96;
  v2(12,3)=(-((T15*T15*(-(y(9)*T194))-(-(y(9)*T111))*(T15*T111+T15*T111))/(T15*T15*T15*T15)));
  v2(13,1)=3;
  v2(13,2)=115;
  v2(13,3)=(T45*T45*(-(params(6)*exp(x(it_, 2))*params(5)*(T130+exp(y(7))*exp(y(7))*getPowerDeriv(exp(y(7)),params(7),2))))-(-(params(6)*exp(x(it_, 2))*params(5)*T130))*(T45*exp(x(it_, 2))*params(5)*T130+T45*exp(x(it_, 2))*params(5)*T130))/(T45*T45*T45*T45);
  v2(14,1)=3;
  v2(14,2)=150;
  v2(14,3)=(-((-T111)/(T15*T15)));
  v2(15,1)=3;
  v2(15,2)=99;
  v2(15,3)=  v2(14,3);
  v2(16,1)=3;
  v2(16,2)=152;
  v2(16,3)=1;
  v2(17,1)=3;
  v2(17,2)=135;
  v2(17,3)=  v2(16,3);
  v2(18,1)=3;
  v2(18,2)=259;
  v2(18,3)=((-(params(6)*exp(x(it_, 2))*params(5)*T130))*T45*T45-(-(params(6)*T45))*(T45*exp(x(it_, 2))*params(5)*T130+T45*exp(x(it_, 2))*params(5)*T130))/(T45*T45*T45*T45);
  v2(19,1)=3;
  v2(19,2)=123;
  v2(19,3)=  v2(18,3);
  v2(20,1)=3;
  v2(20,2)=267;
  v2(20,3)=(T45*T45*(-(params(6)*T45))-(-(params(6)*T45))*(T45*T45+T45*T45))/(T45*T45*T45*T45);
  v2(21,1)=3;
  v2(21,2)=296;
  v2(21,3)=(-(exp(x(it_, 4))*params(4)));
  v2(22,1)=3;
  v2(22,2)=143;
  v2(22,3)=  v2(21,3);
  v2(23,1)=3;
  v2(23,2)=305;
  v2(23,3)=(-(y(8)*exp(x(it_, 4))*params(4)));
  v2(24,1)=4;
  v2(24,2)=96;
  v2(24,3)=y(8)*(-exp(y(6)));
  v2(25,1)=4;
  v2(25,2)=115;
  v2(25,3)=y(8)*(-T65);
  v2(26,1)=4;
  v2(26,2)=132;
  v2(26,3)=(-exp(y(6)));
  v2(27,1)=4;
  v2(27,2)=98;
  v2(27,3)=  v2(26,3);
  v2(28,1)=4;
  v2(28,2)=133;
  v2(28,3)=(-T65);
  v2(29,1)=4;
  v2(29,2)=116;
  v2(29,3)=  v2(28,3);
  v2(30,1)=4;
  v2(30,2)=170;
  v2(30,3)=exp(y(10));
  v2(31,1)=4;
  v2(31,2)=136;
  v2(31,3)=  v2(30,3);
  v2(32,1)=4;
  v2(32,2)=172;
  v2(32,3)=y(8)*exp(y(10));
  v2(33,1)=4;
  v2(33,2)=259;
  v2(33,3)=y(8)*(-T65);
  v2(34,1)=4;
  v2(34,2)=123;
  v2(34,3)=  v2(33,3);
  v2(35,1)=4;
  v2(35,2)=260;
  v2(35,3)=(-T65);
  v2(36,1)=4;
  v2(36,2)=141;
  v2(36,3)=  v2(35,3);
  v2(37,1)=4;
  v2(37,2)=267;
  v2(37,3)=y(8)*(-T65);
  v2(38,1)=4;
  v2(38,2)=277;
  v2(38,3)=T177;
  v2(39,1)=4;
  v2(39,2)=124;
  v2(39,3)=  v2(38,3);
  v2(40,1)=4;
  v2(40,2)=278;
  v2(40,3)=T176;
  v2(41,1)=4;
  v2(41,2)=142;
  v2(41,3)=  v2(40,3);
  v2(42,1)=4;
  v2(42,2)=285;
  v2(42,3)=T177;
  v2(43,1)=4;
  v2(43,2)=268;
  v2(43,3)=  v2(42,3);
  v2(44,1)=4;
  v2(44,2)=286;
  v2(44,3)=y(8)*(-(exp(y(7))*params(5)*exp(x(it_, 2))*params(4)*exp(x(it_, 4))*(T163*(-(exp(x(it_, 3))*params(3)*(1-exp(x(it_, 5))*params(2))))-(-(exp(x(it_, 3))*params(3)*(1-exp(x(it_, 5))*params(2))))*(T163+T163))/(T163*T163)));
  v2(45,1)=4;
  v2(45,2)=295;
  v2(45,3)=y(8)*(-T65);
  v2(46,1)=4;
  v2(46,2)=125;
  v2(46,3)=  v2(45,3);
  v2(47,1)=4;
  v2(47,2)=296;
  v2(47,3)=(-T65);
  v2(48,1)=4;
  v2(48,2)=143;
  v2(48,3)=  v2(47,3);
  v2(49,1)=4;
  v2(49,2)=303;
  v2(49,3)=y(8)*(-T65);
  v2(50,1)=4;
  v2(50,2)=269;
  v2(50,3)=  v2(49,3);
  v2(51,1)=4;
  v2(51,2)=304;
  v2(51,3)=T177;
  v2(52,1)=4;
  v2(52,2)=287;
  v2(52,3)=  v2(51,3);
  v2(53,1)=4;
  v2(53,2)=305;
  v2(53,3)=y(8)*(-T65);
  v2(54,1)=4;
  v2(54,2)=313;
  v2(54,3)=T190;
  v2(55,1)=4;
  v2(55,2)=126;
  v2(55,3)=  v2(54,3);
  v2(56,1)=4;
  v2(56,2)=314;
  v2(56,3)=T189;
  v2(57,1)=4;
  v2(57,2)=144;
  v2(57,3)=  v2(56,3);
  v2(58,1)=4;
  v2(58,2)=321;
  v2(58,3)=T190;
  v2(59,1)=4;
  v2(59,2)=270;
  v2(59,3)=  v2(58,3);
  v2(60,1)=4;
  v2(60,2)=322;
  v2(60,3)=y(8)*(-(exp(y(7))*params(5)*exp(x(it_, 2))*params(4)*exp(x(it_, 4))*(-(exp(x(it_, 3))*params(3)*(-(exp(x(it_, 5))*params(2)))))/T163));
  v2(61,1)=4;
  v2(61,2)=288;
  v2(61,3)=  v2(60,3);
  v2(62,1)=4;
  v2(62,2)=323;
  v2(62,3)=T190;
  v2(63,1)=4;
  v2(63,2)=306;
  v2(63,3)=  v2(62,3);
  v2(64,1)=4;
  v2(64,2)=324;
  v2(64,3)=T190;
  v2(65,1)=5;
  v2(65,2)=96;
  v2(65,3)=exp(y(6));
  v2(66,1)=5;
  v2(66,2)=115;
  v2(66,3)=exp(x(it_, 2))*params(5)*exp(y(7));
  v2(67,1)=5;
  v2(67,2)=172;
  v2(67,3)=(-exp(y(10)));
  v2(68,1)=5;
  v2(68,2)=259;
  v2(68,3)=exp(x(it_, 2))*params(5)*exp(y(7));
  v2(69,1)=5;
  v2(69,2)=123;
  v2(69,3)=  v2(68,3);
  v2(70,1)=5;
  v2(70,2)=267;
  v2(70,3)=exp(x(it_, 2))*params(5)*exp(y(7));
  v2(71,1)=6;
  v2(71,2)=1;
  v2(71,3)=(-(params(3)*exp(x(it_-1, 3))*(-exp(y(1)))));
  v2(72,1)=6;
  v2(72,2)=20;
  v2(72,3)=T128;
  v2(73,1)=6;
  v2(73,2)=38;
  v2(73,3)=T145;
  v2(74,1)=6;
  v2(74,2)=21;
  v2(74,3)=  v2(73,3);
  v2(75,1)=6;
  v2(75,2)=58;
  v2(75,3)=(-(params(3)*exp(x(it_-1, 3))*exp(y(4))));
  v2(76,1)=6;
  v2(76,2)=172;
  v2(76,3)=exp(y(10));
  v2(77,1)=6;
  v2(77,2)=191;
  v2(77,3)=(-exp(y(11)));
  v2(78,1)=6;
  v2(78,2)=254;
  v2(78,3)=T128;
  v2(79,1)=6;
  v2(79,2)=33;
  v2(79,3)=  v2(78,3);
  v2(80,1)=6;
  v2(80,2)=255;
  v2(80,3)=T145;
  v2(81,1)=6;
  v2(81,2)=51;
  v2(81,3)=  v2(80,3);
  v2(82,1)=6;
  v2(82,2)=267;
  v2(82,3)=T128;
  v2(83,1)=6;
  v2(83,2)=271;
  v2(83,3)=(-(params(3)*exp(x(it_-1, 3))*(-exp(y(1)))));
  v2(84,1)=6;
  v2(84,2)=16;
  v2(84,3)=  v2(83,3);
  v2(85,1)=6;
  v2(85,2)=272;
  v2(85,3)=T128;
  v2(86,1)=6;
  v2(86,2)=34;
  v2(86,3)=  v2(85,3);
  v2(87,1)=6;
  v2(87,2)=273;
  v2(87,3)=T145;
  v2(88,1)=6;
  v2(88,2)=52;
  v2(88,3)=  v2(87,3);
  v2(89,1)=6;
  v2(89,2)=274;
  v2(89,3)=(-(params(3)*exp(x(it_-1, 3))*exp(y(4))));
  v2(90,1)=6;
  v2(90,2)=70;
  v2(90,3)=  v2(89,3);
  v2(91,1)=6;
  v2(91,2)=285;
  v2(91,3)=T128;
  v2(92,1)=6;
  v2(92,2)=268;
  v2(92,3)=  v2(91,3);
  v2(93,1)=6;
  v2(93,2)=286;
  v2(93,3)=T157;
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),7,324);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],7,5832);
end
end
end
end
