function [r_curve,norm_Q] = BezierCurvature(P,t,degree,options)
%BezierCurvature computes the curvature radius and the normal direction of
%the curve.
% input:
% - P, control points of Bezier
% - t, time vector of Bezier
% - degree, degree of the Bezier curve (should be at least 3 to compute
% curvature
% - options: "I should check if it works as I want with special options..."

switch nargin
    case 2
       options = t;
    case 3
      if isstruct(degree)  
        options = degree;
      else
          options = struct;
      end
end

if ~isfield(options,'t') && nargin<3
    if isfield(options,'npts')
        t = linspace(0,1,options.npts);
    else
        t = linspace(0,1,250);
    end
elseif isfield(options,'t') && nargin<3
    t = options.t;
end
if isfield(options,'degree')
    degree = options.degree;
end


dQ = BezierDerivative(degree, t, P,1);
ddQ = BezierDerivative(degree, t, P,2);

dQr = c2m(dQ);
ddQr = c2m(ddQ);
num = cross([dQr,zeros(length(dQr),1)],[ddQr,zeros(length(dQr),1)]);
num = num(:,3);
denom = vecnorm(dQr,2,2).^3;

r_curve = denom./num;
norm_Q = 1i*dQ;
end