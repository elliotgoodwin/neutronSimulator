%%%%%%%%%%%%%%%%%% Variation with thickness of shielding %%%%%%%%%%%%%%%%%%
% A script that uses monte-carlo methods to track the path of thermal 
% neutrons through shielding and compares the percentages absorbed,
% transmitted and reflected for different thickness of three materials
% NB c.g.s units used throughout the simulation
% Elliot Goodwin - 9621958

clear all;  close all;  format long;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Avagadro's number
N_A = 6.0221e23;
    
% water parameters
D_w = 1;
A_w = 0.6652e-24;
S_w = 103.0e-24;
M_w = 18.0153;

% lead parameters
D_l = 11.35;
A_l = 0.158e-24;
S_l = 11.221e-24;
M_l = 207.2;

% graphite parameters
D_g = 1.67;
A_g = 0.0045e-24;
S_g = 4.74e-24;
M_g = 12.011;

% mean free path and absorption probabilities
n_w = D_w*N_A/M_w;
AX_w = A_w*n_w;
SX_w = S_w*n_w;
TX_w = AX_w + SX_w;

absorbprob_w = AX_w/TX_w;
mfp_w = 1./TX_w;

n_l = D_l*N_A/M_l;
AX_l = A_w*n_l;
SX_l = S_l*n_l;
TX_l = AX_l + SX_l;

absorbprob_l = AX_l/TX_l;
mfp_l = 1./TX_l;

n_g = D_g*N_A/M_g;
AX_g = A_g*n_g;
SX_g = S_g*n_g;
TX_g = AX_g + SX_g;

absorbprob_g = AX_g/TX_g;
mfp_g = 1./TX_g;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% to find how the number of particles that undergo each process varies with
% thickness of the shielding, the simulation must be repeated for the three
% different shielding materials, with different values of thickness

% number of neutrons simulated for each thickness
k = 5000;

% run simulation ten times 10 times to allow a standard deviation and hence
% uncertainty in the number of neutrons that underwent each process to be
% found

% choose an appropriate interval for independent variable (thickness of
% shield)
interval = 0.5;

for n = 1:10
    
    % run the simulation for differnt shielding thickness from 0 - 20 cm in
    % intervals defined above
    for m = 1:(20./interval)
        
        % save the number of neutrons that underwent each process to an
        % array titled ij where
        % i = {a, r, t}
        % j = {w, l, g}
        % a - absorbed
        % r - reflected
        % t - transmitted
        % w - water
        % l - lead
        % g - graphite
        
        [aw(m, n), rw(m, n), tw(m, n)] = counts(k, mfp_w, absorbprob_w, interval.*m);
        [al(m, n), rl(m, n), tl(m, n)] = counts(k, mfp_l, absorbprob_l, interval.*m);
        [ag(m, n), rg(m, n), tg(m, n)] = counts(k, mfp_g, absorbprob_g, interval.*m);

    end

end

% create vectors containing data to be plotted
% x - path length
% y - percentage of neutrons that have undergone a particular process
% e - error in the percentage of neutrons that have undergone a particular
% process

x = [interval:interval:20]';

% water
Aw = mean(aw, 2);
% error in number of particles
e_Aw = std(aw, 0, 2);
% percentage of particles
Aw_frac = Aw./k';
% error in percentage of particles
e_Aw = e_Aw./k;


Rw = mean(rw, 2);
e_Rw = std(rw, 0, 2);
e_Rw = e_Rw./k;
Rw_frac = Rw./k';

Tw = mean(tw, 2);
Tw_frac = Tw./k';
e_Tw = std(tw, 0, 2);
e_Tw = e_Tw./k;


% lead
Al = mean(al, 2);
e_Al = std(al, 0, 2);
e_Al = e_Al./k;
Al_frac = Al./k';

Rl = mean(rl, 2);
e_Rl = std(rl, 0, 2);
e_Rl = e_Rl./k;
Rl_frac = Rl./k';

Tl = mean(tl, 2);
e_Tl = std(tl, 0, 2);
e_Tl = e_Tl./k;
Tl_frac = Tl./k';

% graphite
Ag = mean(ag, 2);
e_Ag = std(ag, 0, 2);
e_Ag = e_Ag./k;
Ag_frac = Ag./k';

Rg = mean(rg, 2);
e_Rg = std(rg, 0, 2);
e_Rg = e_Rg./k;
Rg_frac = Rg./k';

Tg = mean(tg, 2);
e_Tg = std(tg, 0, 2);
e_Tg = e_Tg./k;
Tg_frac = Tg./k';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Write data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% append results to a vector to be written out to a .txt file
water_data = [x, Aw, Aw_frac, e_Aw, Tw, Tw_frac, e_Tw, Rw, Rw_frac, e_Rw];

% open a file for data to be written to
textdata = fopen('thickness_variation_water.txt','w');
% check file is open
if textdata < 0 
    error('Cannot write to file.');
else
    
    % if file is open, write data to .txt file
    fprintf(textdata, 'Neutron transmission through shielding of varying thickness:\r\n\r\n');
    fprintf(textdata, 'Material:  Water \r\n');
    fprintf(textdata, 'Total number of neutrons: %4.0f\r\n\r\n', [k]);
    fprintf(textdata, '%9s %26s %34s %32s\r\n\r\n', '', 'Absorbed',...
        'Transmitted', 'Reflected');
    fprintf(textdata, '%15s %8s %12s %8s %11s %12s %8s %11s %12s %8s\r\n',...
        'Thickness  [cm]',...
        'Number', 'Percentage', 'Error',...
        'Number', 'Percentage', 'Error',...
        'Number', 'Percentage', 'Error');
    fprintf(textdata, '%9.2f %13.2f %11.6f %13.6f %6.2f %12.6f %13.6f %6.2f %12.6f %13.6f\r\n', water_data');
    fclose(textdata);
    
end


% append results to a vector to be written out to a .txt file
lead_data = [x, Al, Al_frac, e_Al, Tl, Tl_frac, e_Tl, Rl, Rl_frac, e_Rl];

% open a file for data to be written to
textdata = fopen('thickness_variation_lead.txt','w');
% check file is open
if textdata < 0 
    error('Cannot write to file.');
else
    
    % if file is open, write data to .txt file
    fprintf(textdata, 'Neutron transmission through shielding of varying thickness:\r\n\r\n');
    fprintf(textdata, 'Material:  Lead\r\n');
    fprintf(textdata, 'Total number of neutrons: %4.0f\r\n\r\n', [k]);
    fprintf(textdata, '%9s %26s %34s %32s\r\n\r\n', '', 'Absorbed',...
        'Transmitted', 'Reflected');
    fprintf(textdata, '%15s %8s %12s %8s %11s %12s %8s %11s %12s %8s\r\n', 'Thickness  [cm]',...
        'Number', 'Percentage', 'Error',...
        'Number', 'Percentage', 'Error',...
        'Number', 'Percentage', 'Error');
    fprintf(textdata, '%9.2f %13.2f %11.6f %13.6f %6.2f %12.6f %13.6f %6.2f %12.6f %13.6f\r\n', lead_data');
    fclose(textdata);
    
end


% append results to a vector to be written out to a .txt file
graphite_data = [x, Ag, Ag_frac, e_Ag, Tg, Tg_frac, e_Tg, Rg, Rg_frac, e_Rg];

% open a file for data to be written to
textdata = fopen('thickness_variation_graphite.txt','w');
% check file is open
if textdata < 0 
    error('Cannot write to file.');
else
    
    % if file is open, write data to .txt file
    fprintf(textdata, 'Neutron transmission through shielding of varying thickness:\r\n\r\n');
    fprintf(textdata, 'Material:  Water \r\n');
    fprintf(textdata, 'Total number of neutrons: %4.0f\r\n\r\n', [k]);
    fprintf(textdata, '%9s %26s %34s %32s\r\n\r\n', '', 'Absorbed',...
        'Transmitted', 'Reflected');
    fprintf(textdata, '%15s %8s %12s %8s %11s %12s %8s %11s %12s %8s\r\n', 'Thickness  [cm]',...
        'Number', 'Percentage', 'Error',...
        'Number', 'Percentage', 'Error',...
        'Number', 'Percentage', 'Error');
    fprintf(textdata, '%9.2f %13.2f %11.6f %13.6f %6.2f %12.6f %13.6f %6.2f %12.6f %13.6f\r\n', graphite_data');
    fclose(textdata);
    
end



