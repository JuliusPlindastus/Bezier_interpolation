clear
close all
clc

addpath(genpath(pwd));


%% Quick check of Alex's data
load("ms005mpt001_f.mat")

i_time = 120;

QXY = kinematics(i_time).yRawmm;

colors = parula(length(QXY));

figure
scatter(QXY(:, 1), QXY(:, 2), 36, colors, 'filled');
QXY = complex(QXY(:, 1), QXY(:, 2));
%% Bezier fit
options.degree = 4; %By default: 3.
options.anchor_start = true;
options.anchor_last = true;
options.optim = false;%false, true or 'tangent';
%
[Qz,Pz] = BezierFit(QXY.',options);

figure
hold on
plot(QXY,'.r')
plot(Qz,'g','LineWidth',2);
plot(Pz,'-k*');
xlabel('x');ylabel('y');
grid on;
axis equal


[r_c,norm_Q] = BezierCurvature(Pz,options);
Theta_normQ = angle(norm_Q);
r_c(abs(r_c)<10^-5) = 10^-5;
norm = exp(1i.*Theta_normQ)./r_c;
norm(norm > 10) = 10;


quiver(real(Qz),imag(Qz),real(norm),imag(norm),'ShowArrowHead','off','Color',rgb('MediumTurquoise'))

%Just a small example how to reconstruct the curve from the control points
%Pz
options.npts=1000;
Pz2 = Pz+[0;-1i*imag(Pz(2));-20-30i;-50;0];
Qz2 = BezierConstruction(Pz2,options);

figure
hold on
plot(QXY,'.r')
plot(Qz2,'g','LineWidth',2);
plot(Pz,'--k');
plot(Pz2,'-k*');
xlabel('x');ylabel('y');
grid on;
axis equal
