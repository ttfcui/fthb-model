function [residual, g1, g2, g3] = durswitching_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 7, 1);

%
% Model equations
%

T11 = exp(y(1))^params(7);
T12 = 1/T11;
T41 = exp(x(2))*params(5)*exp(y(2))^params(7);
T61 = exp(y(2))*params(5)*exp(x(2))*params(4)*exp(x(4))*(1-exp(x(5))*params(2))/(exp(x(3))*params(3));
lhs =T12;
rhs =(T12-y(3))/(params(1)*exp(x(3))*params(3));
residual(1)= lhs-rhs;
lhs =y(4);
rhs =(exp(x(3))*params(3)-1+exp(x(5))*params(2))/(exp(x(3))*params(3));
residual(2)= lhs-rhs;
lhs =params(6)/T41;
rhs =y(4)/T11+y(3)*(exp(x(4))*params(4)-y(4));
residual(3)= lhs-rhs;
residual(4) = y(3)*(exp(y(5))-exp(y(1))-T61);
lhs =y(7);
rhs =exp(y(5))-exp(y(1))-exp(x(2))*params(5)*exp(y(2));
residual(5)= lhs-rhs;
lhs =exp(y(5));
rhs =exp(x(3))*params(3)*(exp(y(5))-exp(y(1)))+exp(y(2))*params(5)*exp(x(2))*params(3)*exp(x(3))*y(4)+exp(y(6));
residual(6)= lhs-rhs;
lhs =y(6);
rhs =(1-params(8))*params(9)+y(6)*params(8)+x(1);
residual(7)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(7, 7);

  %
  % Jacobian matrix
  %

T90 = exp(y(1))*getPowerDeriv(exp(y(1)),params(7),1);
T93 = (-T90)/(T11*T11);
  g1(1,1)=T93-T93/(params(1)*exp(x(3))*params(3));
  g1(1,3)=(-((-1)/(params(1)*exp(x(3))*params(3))));
  g1(2,4)=1;
  g1(3,1)=(-((-(y(4)*T90))/(T11*T11)));
  g1(3,2)=(-(params(6)*exp(x(2))*params(5)*exp(y(2))*getPowerDeriv(exp(y(2)),params(7),1)))/(T41*T41);
  g1(3,3)=(-(exp(x(4))*params(4)-y(4)));
  g1(3,4)=(-(T12-y(3)));
  g1(4,1)=y(3)*(-exp(y(1)));
  g1(4,2)=y(3)*(-T61);
  g1(4,3)=exp(y(5))-exp(y(1))-T61;
  g1(4,5)=y(3)*exp(y(5));
  g1(5,1)=exp(y(1));
  g1(5,2)=exp(x(2))*params(5)*exp(y(2));
  g1(5,5)=(-exp(y(5)));
  g1(5,7)=1;
  g1(6,1)=(-(exp(x(3))*params(3)*(-exp(y(1)))));
  g1(6,2)=(-(exp(y(2))*params(5)*exp(x(2))*params(3)*exp(x(3))*y(4)));
  g1(6,4)=(-(exp(y(2))*params(5)*exp(x(3))*params(3)*exp(x(2))));
  g1(6,5)=exp(y(5))-exp(x(3))*params(3)*exp(y(5));
  g1(6,6)=(-exp(y(6)));
  g1(7,6)=1-params(8);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],7,49);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],7,343);
end
end
end
end
