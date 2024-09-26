clear
close all
clc

addpath(genpath(pwd));


%% Quick check of Gaetan's data
load("raw_event_with_reconstruct_ms_209_mpt040.mat")
target_value = 2.5;
tolerance = 0.5*5e-4;

time_condition = abs(t_evs - target_value) < tolerance;
XY = [x_evs, y_evs];
XY_tointerpolate = XY(time_condition);

% figure
% plot(XY_tointerpolate, 'x')
% hold on;

field_names = fieldnames(defraw{1});
num_entries = length(defraw);

data = cell(num_entries, length(field_names));

for i = 1:length(field_names)
    field_data = arrayfun(@(s) s.(field_names{i}), [defraw{:}], 'UniformOutput', false);
    data(:, i) = field_data;
end

% Convert the cell array to a table
defraw_table = cell2table(data, 'VariableNames', field_names);

% Convert columns with single numeric values to double arrays
for i = 1:length(field_names)
    if size(defraw_table.(field_names{i}),2)>1
        defraw_table.(field_names{i}) =  num2cell(defraw_table.(field_names{i}), 2);
    end
end

i_time = find(abs(defraw_table.t  - target_value) < tolerance);



% Plot the retrieved x values
% hold on% figure
% plot(defraw_table{i_time,"x"}{:},defraw_table{i_time,"y"}{:}, 'o');
% defraw_table.x{defraw_table.t==2.5}

%% Bezier fit
QXY = complex(defraw_table{i_time,"x"}{:},defraw_table{i_time,"y"}{:});
options.degree = 4; %By default: 3.
options.anchor_start = true;
options.anchor_last = true;
options.optim = false;%false, true or 'tangent';
%
[Qz,Pz] = BezierFit(QXY,options);

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

options.npts=1000;
Qz2 = BezierConstruction(Pz,options);
