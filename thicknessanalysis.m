% Thickness variation analysis
% a script that reads in data containing data from a .txt file to analyse
% how the scattering, absorption and reflection of thermal neutrons
% incident upon shielding varies with thickness

clear all; close all; format long;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read in the percentage and errors of neutrons that underwent each process
% from .txt files
temp = importdata('thickness_variation_water.txt', ' ', 8);
data = temp.data;
xw = data(:, 1);
Aw = data(:, 3);
e_Aw = data(:, 4);
Tw = data(:, 6);
e_Tw = data(:, 7);
Rw = data(:, 9);
e_Rw = data(:, 10);
clear temp data;

temp = importdata('thickness_variation_lead.txt', ' ', 8);
data = temp.data;
xl = data(:, 1);
Al = data(:, 3);
e_Al = data(:, 4);
Tl = data(:, 6);
e_Tl = data(:, 7);
Rl = data(:, 9);
e_Rl = data(:, 10);
clear temp data;

temp = importdata('thickness_variation_graphite.txt', ' ', 8);
data = temp.data;
xg = data(:, 1);
Ag = data(:, 3);
e_Ag = data(:, 4);
Tg = data(:, 6);
e_Tg = data(:, 7);
Rg = data(:, 9);
e_Rg = data(:, 10);
clear temp data;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plot graphs showing the percentages of neutrons that underwent each
% process, for the three shielding materials
figure(1);
errorbar(xw, Aw, e_Aw, 'b.'); hold on;
errorbar(xw, Rw, e_Rw, 'g.');   hold on;
errorbar(xw, Tw, e_Tw, 'r.');    hold on;
grid minor;   hold on;
legend('Aborbed', 'Reflected', 'Transmitted', 'Location', 'Best'); 
xlabel('Shield thickness [cm]'); hold on;
ylabel('Percentage');   hold off;


figure(2);
errorbar(xl, Al, e_Al, 'b.'); hold on;
errorbar(xl, Rl, e_Rl, 'g.');   hold on;
errorbar(xl, Tl, e_Tl, 'r.');    hold on;
grid minor;   hold on;
legend('Aborbed', 'Reflected', 'Transmitted', 'Location', 'Best'); 
xlabel('Shield thickness [cm]'); hold on;
ylabel('Percentage');   hold off;


figure(3);
errorbar(xg, Ag, e_Ag, 'b.'); hold on;
errorbar(xg, Rg, e_Rg, 'g.');   hold on;
errorbar(xg, Tg, e_Tg, 'r.');    hold on;
grid minor;   hold on;
legend('Aborbed', 'Reflected', 'Transmitted', 'Location', 'Best'); 
xlabel('Shield thickness [cm]'); hold on;
ylabel('Percentage');   hold off;

% save figures
saveas(figure(1), 'water_histories', 'png');
saveas(figure(2), 'lead_histories', 'png');
saveas(figure(3), 'graphite_histories', 'png');

% perform linear fits to the log of percentage absorbed and use the
% gradient to find an estimate for the characteristic attenuation length
% for each of the three shielding materials

% water
y = log(Tw);
% error in log(% of particles)
e = e_Tw./Tw;

% remove infinities
Z = zeros(size(xw, 1), 3);
Z(:, 1) = xw(:, 1);
Z(:, 2) = y(:, 1);
Z(:, 3) = e(:, 1);

Z(any(isinf(Z), 2), :) = [];
x = Z(:, 1);
y = Z(:, 2);
e = Z(:, 3);

[p, s] = polyfit(x, y, 1);
f = p(1).*x + p(2);

figure(4);
errorbar(x, y, e, 'o');    hold on;
plot(x, f); hold on;
grid minor; hold on;
xlabel('Shield thickness  [cm]');   hold on;
ylabel('log(r)');   hold on;
legend('Data', 'Linear fit');   hold off;

% construct covariance matrix to extract error on fitting parameters
covm = mean(e).*inv(s.R)*inv(s.R)';
e_gradient = sqrt(covm(1, 1));
e_intercept = sqrt(covm(2, 2));

% estimate mean free path from gradient
mfp = -1./p(1);
maxmfp = (-1./(p(1)-e_gradient));
e_mfp = mfp - maxmfp;
% chi squared analysis
chi2 = ((y - f)./e).^2;
chi2w = sum(chi2)./(s.df);

% summary of results 
mfp_w = [p(1); e_gradient; mfp ; e_mfp; chi2w];



% lead
y = log(Tl);
% error in log(% of particles)
e = e_Tl./Tl;

% remove infinities
Z = zeros(size(xl, 1), 3);
Z(:, 1) = xl(:, 1);
Z(:, 2) = y(:, 1);
Z(:, 3) = e(:, 1);

Z(any(isinf(Z), 2), :) = [];
x = Z(:, 1);
y = Z(:, 2);
e = Z(:, 3);

[p, s] = polyfit(x, y, 1);
f = p(1).*x + p(2);

figure(5);
errorbar(x, y, e, 'o');    hold on;
plot(x, f); hold on;
grid minor; hold on;
xlabel('Shield thickness  [cm]');   hold on;
ylabel('log(r)');   hold on;
legend('Data', 'Linear fit');   hold off;

% construct covariance matrix to extract error on fitting parameters
covm = mean(e)*inv(s.R)*inv(s.R)';
e_gradient = sqrt(covm(1, 1));
e_intercept = sqrt(covm(2, 2));

% estimate mean free path from gradient
mfp = -1./p(1);
maxmfp = (-1./(p(1)-e_gradient));
e_mfp = mfp - maxmfp;
% chi squared analysis
chi2 = ((y - f)./e).^2;
chi2l = sum(chi2)./(s.df);

% summary of data
mfp_l = [p(1); e_gradient; mfp ; e_mfp; chi2l];



% graphite
y = log(Tg);
% error in log(% of particles)
e = e_Tg./Tg;

% remove infinities
Z = zeros(size(xg, 1), 3);
Z(:, 1) = xg(:, 1);
Z(:, 2) = y(:, 1);
Z(:, 3) = e(:, 1);

Z(any(isinf(Z), 2), :) = [];
x = Z(:, 1);
y = Z(:, 2);
e = Z(:, 3);

[p, s] = polyfit(x, y, 1);
f = p(1).*x + p(2);

figure(6);
errorbar(x, y, e, 'o');    hold on;
plot(x, f); hold on;
grid minor; hold on;
xlabel('Shield thickness  [cm]');   hold on;
ylabel('log(r)');   hold on;
legend('Data', 'Linear fit');   hold off;

% construct covariance matrix to extract error on fitting parameters
covm = mean(e)*inv(s.R)*inv(s.R)';
e_gradient = sqrt(covm(1, 1));
e_intercept = sqrt(covm(2, 2));

% estimate mean free path from gradient
mfp = -1./p(1);
maxmfp = (-1./(p(1)-e_gradient));
e_mfp = mfp - maxmfp;
% chi squared analysis
chi2 = ((y - f)./e).^2;
chi2g = sum(chi2)./(s.df);
mfp_g = [p(1); e_gradient; mfp ; e_mfp; chi2g];


% print out summary data to text file

% open a file for data to be written to
textdata = fopen('thickness_variation_summary.txt','w');
%check file is open
if textdata < 0 
    error('Cannot write to file.');
else
    
    % if file is open, write data to .txt file
    fprintf(textdata, 'Neutron transmission through shielding of varying thickness:\r\n\r\n');
    fprintf(textdata, 'Total number of neutrons: %4.0f\r\n\r\n', 5000);
    fprintf(textdata, '%12s %19s %16s %24s %15s %12s\r\n\r\n',...
        'Material', 'Gradient  [/cm]', 'Error  [/cm]',...
        'Mean free path [cm]', 'Error  [cm]', 'Chi-square');
    fprintf(textdata, '%12s %19.4f %16s %24.4f %15.4f %12.4f\r\n', 'Water', mfp_w');
    fprintf(textdata, '%12s %19.4f %16s %24.4f %15.4f %12.4f\r\n', 'Lead', mfp_l');
    fprintf(textdata, '%12s %19.4f %16s %24.4f %15.4f %12.4f\r\n', 'Graphite', mfp_g');
    fclose(textdata);
end

% save figures

saveas(figure(4), 'water_linear', 'png');
saveas(figure(5), 'lead_linear', 'png');
saveas(figure(6), 'graphite_linear', 'png');

