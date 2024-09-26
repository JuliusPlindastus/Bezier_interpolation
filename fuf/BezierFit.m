function [Qzu,Pz] = BezierFit(QXY,options)
% Subject: Interpolate 2D curves with a bunch of known points using Bezier
% functions
% Date:     10 Februrary 2021
% File:     BezierInterpolation.m
% By:       Dr. Julien Deparday
% Source:   Dr. Deparday's thesis and also https://pomax.github.io/bezierinfo/#curvefitting
%
% Input:    QXY:    position of known points with x and y components in columns or rows or as complex arrays
%           options:        structure which contains all additional details
% tuning our plots.
%           options.degree: degree of Bezier polynomial. By default: 3.
%           - options.npts: number of points used to plot the final curve, Qzu.
%           - options.t: if curvilinear distance between known points also
%           known, input it here. options.t must be in an ascending order
%           between 0 and 1.
%           - options.npts_int: number of points used to find the positions
%           of the known points if options.t is not given. By default: 1000.
%           - options.anchor_start: impose position of the first point.
%           - options.anchor_last: impose position of the last point.
%           For both anchors, I don't know how it will behave if t is
%           different from 0 or 1 for the last ones if we force the
%           position of these points...
%           - options.optim:  It can be false, true or 'tangent' to force specific boundary conditions. true by default if options.t not entered
%           - options.optimtang: Force the tangent of the end points to
%           be close to tangent between the two last points. (Useful when
%           many points, and we know how the curve would look like)
%
% Output:   Qzu: Interpolated curve
%           Pz: Position of the control points
%
%
% TODO: Do it for ND-curves
%% Check input

if ~exist('options','var')
    options = struct;
end

if ~isfield(options,'degree')
    options.degree = 3; %By default degree 3
end

if ~isfield(options,'npts_int')
    options.npts_int = 10000; %Use to find best position of known points
end

if ~isfield(options,'npts')
    options.npts = 250; %Number of points for output curve
end

% if ~isfield(options,'optim')
%     %     if isfield(options,'t')
%     options.optim = false;
%     %     else
%     %         options.optim = false;
%     %     end
% end

if isfield(options,'t') % See one day to use try/catch function.
    if ~issorted(options.t)
        error('The curvilinear axis, called here t, must be in an ascending order between 0 and 1!')
    end
    if length(options.t)~=length(QXY)
        error('Each value of the curvilinear axis, t, must correspond to a known point. It seems there isn''t the same number of points...')
    end

    if options.t(end) ~= 1
        warning('Ideally the curvilinear axis, called here t, should be between 0 and 1. I''ll do it for you')
        options.t = options.t/options.t(end);
    end
end

if ~isfield(options,'anchor_start')
    options.anchor_start = false;
end
if ~isfield(options,'anchor_last')
    options.anchor_last = false;
end

if ~isfield(options,'optim')
    options.optimtang = false;
end

if ~isfield(options,'optimtang')
    options.optimtang = false;
end

if ~isfield(options,'step')
    options.step = 10^-6;
end
if ~isfield(options,'iter')
    options.iter = 1000;
end
if ~isfield(options,'thres')
    options.thres = 10^-8;
end
if ~isfield(options,'weight')
    options.weight = ones(size(QXY)).';
end


%Check type of matrix for SXY
if isreal(QXY)
    flagreal = true;

    if size(QXY,1)>size(QXY,2)
        Qx = QXY(:,1)';
        Qy = QXY(:,2)';
        Qz = complex(Qx',Qy');
    else
        Qx = QXY(1,:);
        Qy = QXY(2,:);
        Qz = complex(Qx',Qy');
    end
else
    flagreal = false;
    Qx = real(QXY);
    Qy = imag(QXY);
    Qz = QXY.';
end


%% Find "position" on the curve if not given in options
if ~isfield(options,'t') %&& ~options.optim

    dQz = [0; cumsum(abs(diff(Qz)))];
    options.t = dQz'./dQz(end);
end

%% Computing Control points (Pz)
% Optim tangent
if strcmp(options.optim,'tangent')

    tio = options.t;
    Fdp = @(t) optim_tang(t,options,Qz);
    lb=zeros(size(tio));
    ub=ones(size(tio));
    options_optimtang = optimoptions(@lsqnonlin,'Algorithm','Levenberg-Marquardt',...
        'Display','off','FiniteDifferenceType','central','FunValCheck','on',...
        'StepTolerance',1e-10,'OptimalityTolerance',1e-10); %

    topt = lsqnonlin(Fdp,tio,lb,ub,options_optimtang);

    Pzc = Localoptim(topt,options,Qz);

elseif options.optim%if ~strcmp(options.optim,'tangent')
    tio = linspace(0,1,length(QXY));
    Fdp = @(t) optim_time(t,options,Qz);
    lb=zeros(size(tio));
    ub=ones(size(tio));
    options_optim = optimoptions(@lsqnonlin,'Algorithm','Levenberg-Marquardt',...
        'Display','off','FiniteDifferenceType','central','FunValCheck','on',...
        'StepTolerance',1e-10,'OptimalityTolerance',1e-10); %

    topt = lsqnonlin(Fdp,tio,lb,ub,options_optim);

    Pzc = Localoptim(topt,options,Qz);
else
    Pzc = Localoptim(options.t,options,Qz);
end

%% BEZIERS SURFACE CREATION from control points calculated before
% u for Bezier transformation using control points calculated
% just before. u2 can be different from u.
t_reconstruct = linspace(0,1,options.npts);

[B_u2,U2] = BezierMatrixConstruction(options.degree, t_reconstruct);
% position of the known points Q with the found control points
Qzu = U2*B_u2*Pzc;

if flagreal
    Qzu = [real(Qzu),imag(Qzu)];
    Pz = [real(Pzc),imag(Pzc)];
else
    Pz = Pzc;
end


end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 EXTRA FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B_u,U]=BezierMatrixConstruction(deg_u, u)
% Create the Bezier Matrix with Bernstein polynomials in function of the degree we want
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
% Bernstein polynom
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

% Create the matrix [U^n U^(n-1) ... U 1]
%For deg = 3 [U^3 U^2 U 1]
ku1 = repmat(fliplr(0:deg_u),size(u,2),1); %It gives the

u1 = repmat(u',1,deg_u+1);
U = realpow(u1,ku1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pzc,B_u,ydata_sim] = Localoptim(t,options,Qz)

[B_u,U] = BezierMatrixConstruction(options.degree, t);
B_u_trim = B_u;
Qz_trim = Qz;

if options.anchor_start
    B0 = U*B_u_trim(:,1);
    Q0 = B0*Qz_trim(1);
    Qz_trim = Qz_trim-Q0;
    B_u_trim(:,1) = [];
    options.weight(1)=[];
end
if options.anchor_last
    Bend = U*B_u_trim(:,end);
    Qend = Bend*Qz_trim(end);
    Qz_trim = Qz_trim-Qend;
    B_u_trim(:,end) = [];
    options.weight(end)=[];
end

Pzc_trim = pinv(B_u_trim,1e-15)*pinv(U'*U,1e-15)*U'*Qz_trim;

% Get the inverse function
if options.anchor_start
    Pzc = [Qz(1);Pzc_trim];
else
    Pzc = Pzc_trim;
end

if options.anchor_last
    Pzc = [Pzc;Qz(end)];
end

ydata_sim = U*B_u*Pzc;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DiffNorm = optim_time(t,options,ydata)
[~,~,ydata_sim] = Localoptim(t,options,ydata);
DiffNorm = vecnorm(ydata_sim(2:end-1)-ydata(2:end-1),1,2);
weight = options.weight(2:end-1);
DiffNorm(1:2) = 1.5* DiffNorm(1:2); %small tuning to make it better at the extremities
DiffNorm(end-1:end) = 1.5* DiffNorm(end-1:end);
DiffNorm = DiffNorm.*weight;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DiffTang = optim_tang(t,options,ydata)
[Pzc,~,ydata_sim] = Localoptim(t,options,ydata);

dQz1 = ydata(2)-ydata(1);
dPz1 = Pzc(2)-Pzc(1);
dQzend = ydata(end)-ydata(end-1);
dPzend = Pzc(end)-Pzc(end-1);

condition = [abs(angle(dQz1)-angle(dPz1)), abs(angle(dQzend)-angle(dPzend))];
DiffNorm = vecnorm(ydata_sim(2:end-1)-ydata(2:end-1),1,2);

DiffTang = [DiffNorm.*max(condition);condition.'.*max(DiffNorm)]; %max(condition) is a dirty way to scale it with condition values...
end