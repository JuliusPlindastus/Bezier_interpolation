function [B_u,U]=BezierMatrixConstruction(deg_u, u)
%% Create the Bezier Matrix with Bernstein polynomials in function of the degree we want
% Date:     24 march 2015
% File:     BezierMatrixConstruction.m
% By:       Julien Deparday
% Subject:
% Source:
%
% Input:    deg_u, the degree of the curve or side of surface
%
% Output:   B_u, Bezier-Bernstein polynom
%           U,   matrix of the vector u powered at the degrees deg_u
%           At the end, we will use: Q(u) = U.B_u.P(u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bernstein polynom
% General using Bernstein polynoms
% Bi,n=|n| u^i(1-u)^(n-i) = n!/(i!(n-i)!). u^i.(1-u)^(n-i)
%      |i|
%
%Degree 3 should be: U.B with:
% U = [u^3 u^2 u 1]
% B = [ -1  3 -3  1;...
%        3 -6  3  0;...
%       -3  3  0  0;...
%        1  0  0  0];
%
%
B_u=zeros(deg_u+1,deg_u+1);
ku=1:deg_u+1;
for iu=0:deg_u
    % Coeff n!/(i!(n-i)!)
    Cfact_u = nchoosek(deg_u,iu);
    
    %Binomial coefficient when you develop (1-u)^i, the coefficients follow
    %the Pascal triangle:
    %      1
    %     1 1
    %    1 2 1
    %   1 3 3 1
    %  1 4 6 4 1
    % 1 5 10 10 5 1 etc.
    Binomial_Coeff = diag(rot90(pascal(deg_u-iu+1)))';
    if size(Binomial_Coeff,2) ~= size(B_u(:,iu+1),1)
        Binomial_Coeff(numel(B_u(:,iu+1))) = 0; %Add zeros at the end to fit the matrix B_u
    end
    
    %Signs in function of odd or even
    if mod(deg_u-iu,2)==1
        Signs=(-1).^ku; %(-1)^j
    else
        Signs=-1.*(-1).^ku; %(-1)^(j+1)
    end
    
    B_u(:,iu+1) = (Cfact_u.*Signs.*Binomial_Coeff)';
end

%% Create the matrix [U^n U^(n-1) ... U 1]
%For deg = 3 [U^3 U^2 U 1]
ku1 = repmat(fliplr(0:deg_u),size(u,2),1); %It gives the

u1 = repmat(u',1,deg_u+1);
U = realpow(u1,ku1);
end

