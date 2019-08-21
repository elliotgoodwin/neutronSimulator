% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Project 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script that uses monte-carlo methods to track the path of thermal 
% neutrons through shielding
% NB c.g.s units used throughout the simulation
% Elliot Goodwin - 9621958

clear all;  close all;  format long;

% global variables
% Avagadro's number
N_A = 6.0221e23;

% Parameter definitions

    % D - density
    % A - absorbtion cross section
    % S - scattering cross section
    % M - molar mass
    % n - number of particles available for interaction
    % AX - macroscopic absorption cross section
    % SX - macroscopic scatter cross section
    % TX - total interaction cross section

    
% open user interface to set up simulation
prompt = {'Shield material (Water / w, Lead / l, Graphite / g)',...
    'Shielding thickness  [cm]'...
    'Number of particles in simulation'};
defAns = {'Water', '10', '10000'};
params = inputdlg(prompt, 'Set up simulation', 1, defAns);

% if user presses cancel, exit programme
if isempty(params) == 1
    clear all;
    error('User terminated programme.');
end

% extract parameters
material = params{1, 1};
x = str2double(params{2, 1});
k = str2double(params{3, 1});

% check parameters are sensible
% loop if parameters are invalid

% if the string 'material' doesn't match the desired input, open help
% dialog
if strcmp(material, 'l') == 1 | strcmp(material, 'Water') == 1 |...
   strcmp(material, 'w') == 1 | strcmp(material, 'Graphite') == 1 | ...
   strcmp(material, 'g') == 1 | strcmp(material, 'Lead') == 1

    stringcheck = 1;
else 
    stringcheck = 0; 
end

while k <= 0 | rem(k, 1) ~= 0 | x <= 0 | stringcheck ~= 1;

    % display help window
    uiwait(helpdlg({'Error:', '',...
        'Please ensure parameters are input in the following format',...
        'Number of particles', '',...
        'Any positive integer', '', '', ...
        'Shield material', '',...
        'Water / w', '', 'Lead / l', '', 'Graphite, g'}, 'Error'));
    
    % re-ask for parameter input
    prompt = {'Shield material (Water / w, Lead / l, Graphite / g)',...
        'Shielding thickness  [cm]'...
        'Number of particles in simulation'};
    defAns = {'Water', '10', '10000'};
    params = inputdlg(prompt, 'Set up simulation', 1, defAns);

    % if user presses cancel, exit programme
    if isempty(params) == 1
        error('User terminated programme.');
    end

    % extract parameters
    material = params{1, 1};
    x = str2double(params{2, 1});
    k = str2double(params{3, 1});

end
% if the user selected water shielding
if strcmp(material, 'Water') == 1 | strcmp(material, 'w') == 1
    
    % water parameters
    D = 1;
    A = 0.6652e-24;
    S = 103.0e-24;
    M = 18.0153;

% if the user selected lead shielding
elseif strcmp(material, 'Lead') == 1 | strcmp(material, 'l') == 1

    % lead parameters
    D = 11.35;
    A = 0.158e-24;
    S = 11.221e-24;
    M = 207.2;

% if the user selected graphite shielding
elseif strcmp(material, 'Graphite') == 1 | strcmp(material, 'g') == 1
    
    % graphite parameters
    D = 1.67;
    A = 0.0045e-24;
    S = 4.74e-24;
    M = 12.011;
    
end

% calculate probability of absorption
n = D*N_A/M;
AX = A*n;
SX = S*n;
TX = AX + SX;

absorbprob = AX/TX;

% check that this makes physical sense
scatterprob = SX/TX;
if scatterprob + absorbprob ~= 1
    errordlg('Ensure parameters have been input correctly.');
end
   
% calculate mean free path (accounting for scattering and absorption)
mfp = 1/TX;

% check that this makes physical sense
if mfp > 1/AX
    errordlg('Ensure parameters have been input correctly.');   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% simulate paths of k particles
[absorb_count, reflect_count, transmit_count] = counts(k, mfp, absorbprob, x);
[data] = particlehist(k, mfp, absorbprob, x);

% plot particle path
% extract the most 'interesting' (path with the largest number of entries)
% from the particle path structure
for n = 1:k
    
    if data.x(end, n) ~= 0;
    
        r1 = data.x(:, n);
        r2 = data.y(:, n);
        r3 = data.z(:, n);
    
    end

end

% add axes labels

figure(1)
plot3(r1, r2, r3);  hold on;
grid;   hold on;
xlabel('x  [cm]');  hold on;
ylabel('y  [cm]');  hold on;
zlabel('z  [cm]');  hold off;
saveas(figure(1), 'neutronpath', 'png');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% fraction of particles that undergo each process
transmit = transmit_count./k;
reflect= reflect_count./k;
absorb = absorb_count./k;
% absolute error on number of particles for each process
err_transmit = 1./sqrt(transmit_count);
err_absorb = 1./sqrt(absorb_count);
err_reflect = 1./sqrt(reflect_count);
% error on fraction of particles that underwent each process
err_T = (transmit_count + err_transmit)./k - transmit;
err_R = (reflect_count + err_reflect)./k - reflect;
err_A = (absorb_count + err_absorb)./k - absorb;

% ensure errors make physical sense by removing infinities
if err_transmit == inf
    err_transmit = 0;
end
if err_absorb == inf
    err_absorb = 0;
end
if err_reflect == inf
    err_reflect = 0;
end
if err_T == inf
    err_T = 0;
end
if err_R == inf
    err_R = 00;
end
if err_A == inf
    err_A = 0;
end


% write data out to .txt file
% append results to a vector to be written out to a .txt file
T = [transmit, err_T];
R = [reflect, err_T];
A = [absorb, err_A];

% open a file for data to be written to
textdata = fopen('tempdata.txt','w');
% check file is open
if textdata < 0 
    error('Cannot write to file.');
else
    
    % if file is open, write data to .txt file
    fprintf(textdata, 'Neutron transmission through a fixed thickness of shielding\r\n\r\n');
    fprintf(textdata, 'Material:  %s \r\n', [material]);
    fprintf(textdata, 'Thickness: %4.2f cm \r\n', [x]);
    fprintf(textdata, 'Total number of neutrons: %4.0f\r\n\r\n', [k]);
    fprintf(textdata, 'Number reflected:    %4.0f +/- %.4f \r\n', [reflect_count, err_reflect]);
    fprintf(textdata, 'Number absorbed:     %4.0f +/- %.4f \r\n', [absorb_count, err_absorb]);
    fprintf(textdata, 'Number transmitted:  %4.0f +/- %.4f \r\n\r\n', [transmit_count, err_transmit]);
    fprintf(textdata, 'Percent reflected:   %.6f +/- %.6f \r\n', R);
    fprintf(textdata, 'Percent absorbed:    %.6f +/- %.6f \r\n', A);
    fprintf(textdata, 'Percent transmitted: %.6f +/- %.6f \r\n', T);
    fclose(textdata);
    
end



%%%%%%%%%%%%%%%%%%% Characteristic attenuation length %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% find mean free paths from data
% find maximum path length of each neutron before absoption
xdistance = max(data.x(:, :));

% plot a histogram of sampled values
figure(2);
[n, m] = hist(xdistance, 30);
hist(xdistance, 30);    hold on;
xlabel('Path length [cm]');
ylabel('Frequency');
grid;   hold off;
% remove zeroes
m = m(n~=0);
n = n(n~=0);
% perform a linear fit to binned data
[p, st] = polyfit(m, log(n), 1);
f = p(1).*m + p(2);
% error on number in each bin, n ~ 1/sqrt(n)
err = 1./sqrt(n);
figure(3);
% plot data
errorbar(m, log(n), err, 'o');   hold on;
% and line of best fit
plot(m, f); hold on;
grid;   hold on;
xlabel('Path length [cm]'); hold on;
ylabel('log(n)');    hold on;
legend('Data', 'Linear fit');

% construct covariance matrix
covm = mean(err)*inv(st.R)*inv(st.R)';
% error on gradient
err_gradient = sqrt(covm(1,1));
% max value of gradient
max_gradient = p(1) + err_gradient;

% extract mfp and error from the gradient of the fit
atten_length = -1./p(1);
max_atten_length = -1./max_gradient;
err_atten_length = abs(max_atten_length - atten_length);
disp(atten_length);
disp(err_atten_length);


saveas(figure(2), 'pathhistogram', 'png');
saveas(figure(3), 'linearfit_water', 'png');



%%%%%%%%%%%%%%%%%%%%%%%%% Varying neutron number %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% run simulation for number of neutrons = [500:500:10000]
% simulation uses 10 cm water as shielding medium
for n = 1:20

        [a(n), r(n), t(n)] = counts(n.*500, mfp, absorbprob, 10);
        y(n) = ((1./sqrt(t(n)))./(n.*500));

end


x = [500:500:10000];

% plot graph of error variance with shield thickness
figure(4);
plot(x, y, 'bo');   hold on;
grid minor;
xlabel('Number of neutrons');   hold on;
ylabel('Error in percentage of neutrons transmitted');  hold on;

saveas(figure(4), 'neutron_number', 'png');



