function Qz = BezierConstruction(Pz,options)
%% BEZIERS SURFACE CREATION from control points calculated before
% Input:    Pz: control points
%           options: to get the number of points, the degree etc.
% Output:   Interpolated curve Qz

%% Check input
if ~exist('options','var')
    options = struct;
end

if ~isfield(options,'degree')
    options.degree = 3; %By default degree 3
end

if ~isfield(options,'npts')
    options.npts = 250; %Number of points for output curve
end

%Check type of matrix for SXY
if isreal(Pz)
    flagreal = true;

    if size(Pz,1)>size(Pz,2)
        Pzx = Pz(:,1)';
        Pzy = Pz(:,2)';
        Pz = complex(Pzx',Pzy');
    else
        Pzx = Pz(1,:);
        Pzy = Pz(2,:);
        Pz = complex(Pzx',Pzy');
    end
else
    flagreal=false;
end

%% prepare output
t_reconstruct = linspace(0,1,options.npts);

[B_u2,U2] = BezierMatrixConstruction(options.degree, t_reconstruct);
% position of the known points Q with the found control points
Qz = U2*B_u2*Pz;

if flagreal
    Qzu = [real(Qz),imag(Qz)];
end