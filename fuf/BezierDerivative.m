function [Q_der]=BezierDerivative(deg_u, t, Pc,n)
%% Derivate the Bezier as many times as we want
% Date:     23 June 2021
% File:     BezierDerivative.m
% By:       Dr. Julien Deparday
% Subject:
% Source:
%
% Input:    deg_u, the degree of the curve or side of surface
%           t, for the curvilinear axis
%           Pc, control point
%           n, n-th derivative of the Bezier Qz = U*B_u*Pc;
% Output:   B_u, Bezier-Bernstein polynom
%           U,   matrix of the vector u powered at the degrees deg_u
%           At the end, we will use: Q(u) = U.B_u.P(u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if deg_u <= n
    error('The order of the derivation is higher than the order of the Bezier polynom. It won''t work, boss...')
end

% derivate as many times as needed.
dPc = Pc;
for ider=1:n
    dPc = diff(dPc);
end

[B_u_prime,U_prime] = BezierMatrixConstruction(deg_u-n, t);
Q_der = deg_u.*U_prime*B_u_prime*dPc;

