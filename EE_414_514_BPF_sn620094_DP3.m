%------------------------------------------------------------
%% EE-414-514 sn620094 Design Project 3
%------------------------------------------------------------

clearvars
clc
close all

%----------------------------------------------------------------------
% Units 
%----------------------------------------------------------------------

G = 10^9;
Meg = 10^6;
k = 10^+3;
c = 10^-2;
m = 10^-3;
u = 10^-6;
n = 10^-9;

%----------------------------------------------------------------------
% 
%----------------------------------------------------------------------
IFigure = 0;
NF = 32;
dfreq = 1;
df = 1*Meg;
j = 1*j;
theta_P_f0 = 90;
%----------------------------------------------------------------------
% 
%---------------------------------------------------------------------
fls = 21*G;
fhs = 23.4*G;
f0 = 22*G;
BWf = 1*G;

fL = f0 - BWf/2;
fH = f0 + BWf/2;


delta = ((fH - fL)/f0);

Print_Real_Unit('f0',f0,'Hz')
Print_Real_Unit('Bwf',BWf,'Hz')
Print_Real_Unit('fL',fL,'Hz')
Print_Real_Unit('fH',fH,'Hz')
Print_Real('delta',delta)

%%
%----------------------------------------------------------------------
% ILmin = 0.5; % dB
% S21_dB = -ILmin;
% S21min = 10^((S21_dB)/20);
% Print_Real('S21_Min(IL)',S21min,'W/W'); %Good
% 
% S11_max = sqrt(1-abs(S21min)^2); % W/W
% Print_Real('S11_max (IL)',S11_max,'W/W'); %Good
% S11_dB = 20*log10(S11_max); % IL
% 
% Print_Real('S11_dB (IL)',S11_dB,'W/W'); %Good

RLmin = 26;

S11_dB_RL = -RLmin;

S11_max_RL = 10^(((S11_dB_RL)/20));

S21min_RL = sqrt(1-abs(S11_max_RL)^2);
Print_Real('S21_min (RL)',S21min_RL,'W/W'); %Good
S21_dB_RL = 20*log10(S21min_RL); % IL
Print_Real('S21_min (RL)',S21_dB_RL,'dB'); %Good

Ap_dB = abs(S21_dB_RL);
%Ap_dB = 0.02; % Round down
Print_Real('Ap (RIPPLE)',Ap_dB,'dB'); %Good
%----------------------------------------------------------------------
a = 10.6680*m;
b = 4.3180*m;

uc = physconst('LightSpeed');
fc = uc / (2*a);
lambda_fs0 = uc/f0;
lambda_g0 = lambda_fs0/sqrt(1-(fc/f0)^2);

Print_Real_Unit('fc',fc,'Hz')
Print_Real_Unit('lambda_fs0',lambda_fs0,'m')
Print_Real_Unit('lambda_g0',lambda_g0,'m')
%----------------------------------------------------------------------
no = 376.7303;
ZTE = no/sqrt(1-(fc/f0)^2);
Z0 = (2*b/a)*ZTE; %Use
Print_Real_Unit('ZTE',ZTE,'Ohms')
Print_Real_Unit('Z0',Z0,'Ohms')
%%
%----------------------------------------------------------------------

wp = 2*pi*f0;
Ap = 10^(Ap_dB/10); % W/W
Xp = Ap - 1;
epsilon = sqrt(Xp);
%----------------------------------------------------------------------
% FLS
%----------------------------------------------------------------------

ALS = 20;
ALS_W = 10^(ALS/10);
XLS = ALS_W - 1;
Omega_LS = (1/delta)*((fls/f0)-(f0/fls));
NLS = acosh(sqrt(XLS)/epsilon ) / acosh(abs(Omega_LS));
N_LS = ceil(NLS);
NLS_Elements = ceil(N_LS/2);
Print_Real('NLS',NLS)
%----------------------------------------------------------------------
% FHS
%----------------------------------------------------------------------
AHS = 35;
AHS_W = 10^(AHS/10);
XHS = AHS_W - 1;
Omega_HS = (1/delta)*((fhs/f0)-(f0/fhs));
NHS = acosh(sqrt(XHS)/epsilon ) / acosh(abs(Omega_HS));
Print_Real('NHS',NHS)
% N_HS = ceil(NHS) ;
% NHS_Elements = ceil(N_HS/2);
% NP = 5; % Round up twice in example
% %
%%
%Z0 = 50;
Y0 = 1/Z0;
w0 = 2*pi*f0;
Print_Real_Unit('w0',w0/G,'Grad/s')
g0 = 1.0000;
g1 = 0.7563;
g2 = 1.3049;
g3 = 1.5773;
g4 = 1.3049;
g5 = 0.7563;
g6 = 1.0000;


% 
%----------------------------------------------------------------------
% Branch #1 (Shunt)
%----------------------------------------------------------------------
C1 = (1/(w0*delta))*g1*Y0;
L1 = 1/(w0^2*C1);
Print_Break
Print_Real_Unit('C1',C1,'F')
Print_Real_Unit('L1',L1,'H')

%----------------------------------------------------------------------
% Branch #2 (Series)
%----------------------------------------------------------------------
L2 = (1/(w0*delta)) * g2 * Z0;
C2 = 1/(w0^2*L2);
Print_Real_Unit('L2',L2,'H')
Print_Real_Unit('C2',C2,'F')
%----------------------------------------------------------------------
% Branch #3 (Shunt)
%----------------------------------------------------------------------
C3 = (1/(w0*delta))*g3*Y0;
L3 = 1/(w0^2*C3);
Print_Real_Unit('C3',C3,'F')
Print_Real_Unit('L3',L3,'H')
%----------------------------------------------------------------------
% Branch #4 (Series)
%----------------------------------------------------------------------
L4 = (1/(w0*delta)) * g4 * Z0;
C4 = 1/(w0^2*L4);
Print_Real_Unit('L4',L4,'H')
Print_Real_Unit('C4',C4,'F')
%----------------------------------------------------------------------
% Branch #5 (Shunt)
%----------------------------------------------------------------------
C5 = (1/(w0*delta))*g5*Y0;
L5 = 1/(w0^2*C5);
Print_Real_Unit('C5',C5,'F')
Print_Real_Unit('L5',L5,'H')

f_min = 20.5*G;
f_max = 23.5*G;
freq = f_min : df : f_max;
freq = sort(freq);
freq = freq';
% I_f0_Amp = freq_Amp == f0;
%Z0 = 50;
N_Freq = length(freq);
S_Filter = zeros(N_Freq,2,2);
for kk = 1 : N_Freq
fk = freq(kk);
T0 = eye(2); % ?
T1 = EE414_ABCD_Shunt_C1(C1, fk);
T2 = EE414_ABCD_Shunt_L1(L1, fk);
T3 = EE414_ABCD_Series_L2(L2,fk);
T4 = EE414_ABCD_Series_C2(C2,fk);
T5 = EE414_ABCD_Shunt_C3(C3, fk);
T6 = EE414_ABCD_Shunt_L3(L3, fk);
T7 = EE414_ABCD_Series_L4(L4,fk);
T8 = EE414_ABCD_Series_C4(C4,fk);
T9 = EE414_ABCD_Shunt_C5(C5, fk);
T10 = EE414_ABCD_Shunt_L5(L5, fk);
T = T0*T1*T2*T3*T4*T5*T6*T7*T8*T9*T10;
S_Filter(kk, :, :) = ABCD_to_S(T, [Z0, Z0]);
end

%----------------------------------------------------------------------
%
%----------------------------------------------------------------------

S11_Filter = S_Filter(:,1,1);
S11_Filter_Mag = abs(S11_Filter);
S11_Filter_dB = 20*log10(S11_Filter_Mag);
S21_Filter = S_Filter(:,2,1);
S21_Filter_Mag = abs(S21_Filter);
S21_Filter_dB = 20*log10(S21_Filter_Mag);


%----------------------------------------------------------------------
% Impedance Inverter
%----------------------------------------------------------------------
lambda_fs_fL = uc/fL;
lambda_g_fL = lambda_fs_fL/sqrt(1-(fc/fL)^2);
%
lambda_fs_f0 = uc/f0;
lambda_g_f0 = lambda_fs_f0/sqrt(1-(fc/f0)^2);
%
lambda_fs_fH = uc/fH;
lambda_g_fH = lambda_fs_fH/sqrt(1-(fc/fH)^2);
%
Delta_Lambda = (lambda_g_fL - lambda_g_fH)/lambda_g_f0;
%----------------------------------------------------------------------
K01_Z0 = sqrt(((1/2)*pi) * Delta_Lambda / (g0 * g1));
K12_Z0 = (((1/2)*pi) * Delta_Lambda / (sqrt(g1 * g2)));
K23_Z0 = (((1/2)*pi) * Delta_Lambda / (sqrt(g2 * g3)));
K34_Z0 = (((1/2)*pi) * Delta_Lambda / (sqrt(g3 * g4)));
K45_Z0 = (((1/2)*pi) * Delta_Lambda / (sqrt(g4 * g5)));
K56_Z0 = sqrt(((1/2)*pi) * Delta_Lambda / (g5 * g6));
%----------------------------------------------------------------------
% Run HFSS Simulation to find ag, Use 
%----------------------------------------------------------------------
Print_Break
Print_Real_Unit('lambda_fs',lambda_fs_fL,'m')
Print_Real_Unit('lambda_g_fL',lambda_g_fL,'m')
Print_Real_Unit('lambda_fs_f0',lambda_fs_f0,'m')
Print_Real_Unit('lambda_g_f0',lambda_g_f0,'m')
Print_Real_Unit('lambda_fs_fH',lambda_fs_fH,'m')
Print_Real_Unit('lambda_g_fH',lambda_g_fH,'m')
Print_Real('Delta_Lambda',Delta_Lambda)
Print_Break
Print_Real('K01_Z0',K01_Z0)
Print_Real('K12_Z0',K12_Z0)
Print_Real('K23_Z0',K23_Z0)
Print_Real('K34_Z0',K34_Z0)
Print_Real('K45_Z0',K45_Z0)
Print_Real('K56_Z0',K56_Z0)
%% Define the path to the HFSS exported CSV file
csvFilePath = 'C:\Users\s224lab\Desktop\EE-414-514-Design Project 3\S Parameter Plot pec.csv';

% Read the CSV file into a table
data = readtable(csvFilePath);

% Extract the frequency and S11 data from the table
Frequency1 = data.frequency1;  % Assuming the column name in the CSV is 'Frequency'
S11 = data.S11;              % Assuming the column name in the CSV is 'S11'
S21 = data.S21;
NF = 32;
% Calculate the gain (magnitude of S11)
gainS11 = (S11);
gains21 = (S21);


%
% Define the path to the HFSS exported CSV file
csvFilePath = 'C:\Users\s224lab\Desktop\EE-414-514-Design Project 3\S Parameter copper.csv';

% Read the CSV file into a table
data = readtable(csvFilePath);

% Extract the frequency and S11 data from the table
Frequency2 = data.frequency2;  % Assuming the column name in the CSV is 'Frequency'
S11_C = data.S11_C;              % Assuming the column name in the CSV is 'S11'
S21_C = data.S21_C;
NF = 32;
% Calculate the gain (magnitude of S11)
gainS11_C = (S11_C);
gains21_C = (S21_C);


%----------------------------------------------------------------------

IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S11_Filter_dB, 'r', 'linewidth', 6)
hold on
plot(freq/G, S21_Filter_dB, 'g', 'linewidth', 6);
plot(f0/G,interp1(freq/G,S21_Filter_dB,f0/G),'go','linewidth',10);
plot(fL/G,interp1(freq/G,S21_Filter_dB,fL/G),'go','linewidth',10);
plot(fL/G,interp1(freq/G,S11_Filter_dB,fL/G),'ro','linewidth',10);
plot(fH/G,interp1(freq/G,S21_Filter_dB,fH/G),'go','linewidth',10);
plot(fH/G,interp1(freq/G,S11_Filter_dB,fH/G),'ro','linewidth',10);
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{k1} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
title('MATLAB Ideal Lumped Elements')
legend('|{\itS}_{11}|','|{\itS}_{21}|')
axis([f_min/G, f_max/G, -40, 0])
set(gca, 'xtick', f_min/G : 0.1 : f_max/G);
set(gca, 'linewidth', 2.0)
%

IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S11_Filter_dB, 'g', 'linewidth', 6);
hold on
plot(Frequency1, gainS11,'r', 'linewidth',4)
plot(Frequency2, gainS11_C,'b', 'linewidth',4)
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{11} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
%title('MATLAB Ideal Lumped Elements')
legend('Lumped Elements','HFSS w/o Loss','HFSS w/ Loss')
axis([f_min/G, f_max/G, -40, 0])
set(gca, 'xtick', f_min/G : 0.2 : f_max/G);
set(gca, 'linewidth', 2.0)

%----------------------------------------------------------------------
% Iris Sweep Values
%----------------------------------------------------------------------
f0 = 22*G;
%--------------------------------------------------------------------

% Folder_Main_1 =...
%     'G:\EE 414 Microwave Engineering\Design Project 03';
% % Folder_Main_2 =...
% %     'Filters\Waveguide Filters';
% Folder_Main =...
%     [Folder_Main_1, '\', Folder_Main_2];

%--------------------------------------------------------------------

Folder_HFSS = 'C:\Users\s224lab\Desktop\EE-414-514-Design Project 3';
FileName_HFSS_Iris = ...
    'S Parameter Table 4.csv';
File_Loc_HFSS_Iris = ...
    [Folder_HFSS, '\', FileName_HFSS_Iris];


%--------------------------------------------------------------------
 K_Z0 = [ 0.3998, 0.1217, 0.0843, 0.0843, 0.1217, 0.3998]
%--------------------------------------------------------------------

Data_Iris = ...
    readmatrix(File_Loc_HFSS_Iris);
Data_Iris = transpose(Data_Iris);
a_g_HFSS = Data_Iris(2, :);
a_g_HFSS = a_g_HFSS*m;
S11_Iris_Mag = Data_Iris(3, :);
S11_Iris_Ang = Data_Iris(4, :);
S11_Iris = ...
    Polar_2_Rect(S11_Iris_Mag, S11_Iris_Ang);
S21_Iris_Mag = Data_Iris(5, :);
S21_Iris_Ang = Data_Iris(6, :);
S21_Iris = ...
    Polar_2_Rect(S21_Iris_Mag, S21_Iris_Ang);

%--------------------------------------------------------------------

%--------------------------------------------------------------------

Xp_Z0_HFSS = zeros(size(a_g_HFSS));
Xs_Z0_HFSS = zeros(size(a_g_HFSS));
phi_HFSS = zeros(size(a_g_HFSS));
psi_HFSS = zeros(size(a_g_HFSS));
K_Z0_HFSS = zeros(size(a_g_HFSS));
for kk = 1 : length(a_g_HFSS)
    S_HFSS_Iris = [S11_Iris(kk), S21_Iris(kk);
        S21_Iris(kk), S11_Iris(kk)];
    [Xp_Z0_HFSS(kk), Xs_Z0_HFSS(kk), K_Z0_HFSS(kk),...
        phi_HFSS(kk), psi_HFSS(kk)] = ...
        EE414_T_Model_Param(S_HFSS_Iris);
end

%--------------------------------------------------------------------


%--------------------------------------------------------------------

delta_a_g = (a_g_HFSS(end)-a_g_HFSS(1))/100;
a_g_HFSS_2 = a_g_HFSS(1) : delta_a_g : a_g_HFSS(end);
P_K_Z0_vs_ag = spline(a_g_HFSS, K_Z0_HFSS);
K_Z0_HFSS_2 = ppval(P_K_Z0_vs_ag, a_g_HFSS_2);
P_Xs_Z0_vs_ag = spline(a_g_HFSS, Xs_Z0_HFSS);
Xs_Z0_2 = ppval(P_Xs_Z0_vs_ag, a_g_HFSS_2);
P_Xp_Z0_vs_ag = spline(a_g_HFSS, Xp_Z0_HFSS);
Xp_Z0_2 = ppval(P_Xp_Z0_vs_ag, a_g_HFSS_2);
P_phi_vs_ag = spline(a_g_HFSS, phi_HFSS);
phi_2 = ppval(P_phi_vs_ag, a_g_HFSS_2);
P_psi_vs_ag = spline(a_g_HFSS, psi_HFSS);
psi_2 = ppval(P_psi_vs_ag, a_g_HFSS_2);
a_g_HFSS_Plot = a_g_HFSS(1) : delta_a_g : a_g_HFSS(end);
K_Z0_HFSS_Plot = ppval(P_K_Z0_vs_ag, a_g_HFSS_Plot);
phi_Plot = ppval(P_phi_vs_ag, a_g_HFSS_Plot);
%--------------------------------------------------------------------
IFigure = 0;
NF = 32;
%--------------------------------------------------------------------
 
IFigure = IFigure + 1;
figure_max(IFigure)
plot(a_g_HFSS/m, K_Z0_HFSS, 'kx', 'linewidth', 15);
hold on
plot(a_g_HFSS_2/m, K_Z0_HFSS_2, 'r', 'linewidth', 9);
plot(a_g_HFSS/m, K_Z0_HFSS, 'kx', 'linewidth', 15);
hold off
title('Cubic Spline Interpolation')
xlabel('{\ita_g}  (mm)')
ylabel(' {\itK} / {\itZ}_{0}   (\Omega/\Omega)',...
    'VerticalAlignment', 'bottom')
legend(...
    ' HFSS Simulation',...
    ' CS Interpolation',...
    'Location', 'Best')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
grid on
grid minor
set(gca, 'linewidth', 2)
 
%--------------------------------------------------------------------
ag_vs_P_K_Z0 = spline(K_Z0_HFSS,a_g_HFSS);
 a_g = ppval(ag_vs_P_K_Z0, K_Z0);

 
Xs_Z0 = ppval(P_Xs_Z0_vs_ag, a_g);
 Xp_Z0 = ppval(P_Xp_Z0_vs_ag, a_g);
 % psi = atand(Xs_Z0);
 % phi = (-1)*( atand(2*Xp_Z0+Xs_Z0) + psi );
 psi = ppval(P_psi_vs_ag, a_g);
 phi = ppval(P_phi_vs_ag, a_g);
 % K_Z0 = abs(tand(0.5*phi + psi));
 K_Z0 = ppval(P_K_Z0_vs_ag, a_g);
 table(a_g'/m, K_Z0')
a_g_m = a_g/m;
%

theta_k_1 = 180 + (1/2) * (phi(1) + phi(2));
theta_k_2 = 180 + (1/2) * (phi(2) + phi(3));
theta_k_3 = 180 + (1/2) * (phi(3) + phi(4));
theta_k_4 = 180 + (1/2) * (phi(4) + phi(5));
theta_k_5 = 180 + (1/2) * (phi(5) + phi(6));

%
a = 10.6680*m;
b = 4.3180*m;

uc = physconst('LightSpeed');
fc = uc / (2*a);
lambda_fs0 = uc/f0;
lambda_g0 = lambda_fs0/sqrt(1-(fc/f0)^2);

Print_Real_Unit('fc',fc,'Hz')
Print_Real_Unit('lambda_fs0',lambda_fs0,'m')
Print_Real_Unit('lambda_g0',lambda_g0,'m')

% radians

dk_1_lambda_g0 = (theta_k_1/360) ;
dk_2_lambda_g0 = (theta_k_2/360) ;
dk_3_lambda_g0 = (theta_k_3/360) ;
dk_4_lambda_g0 = (theta_k_4/360) ;
dk_5_lambda_g0 = (theta_k_5/360) ;

% dk length in mm / need to find d0 and d6
% d0 = (1/2)*a
% d6 = d0
dk0 = (1/2)*a;
dk1 = dk_1_lambda_g0 * lambda_g0;
dk2 = dk_2_lambda_g0 * lambda_g0;
dk3 = dk_3_lambda_g0 * lambda_g0;
dk4 = dk_4_lambda_g0 * lambda_g0;
dk5 = dk_5_lambda_g0 * lambda_g0;
dk6 = (1/2)*a;

% z  length which is zk in HFSS
t = 1*m;

z1 = 0 + dk0;
z2 = z1 + t + dk1;
z3 = z2 + t +dk2;
z4 = z3 + t + dk3;
z5 = z4 + t + dk4;
z6 = z5 + t + dk5;
z7 = z6 + t + dk6;

L = z7;
dp = dk0;

Print_Real_Unit('phi',phi,'degrees')
Print_Real_Unit('dk0',dk0,'m')

Print_Real_Unit('theta_k_1',theta_k_1,'degrees')
Print_Real('dk_1_lambda_g0',dk_1_lambda_g0)
Print_Real_Unit('dk1',dk1,'m')

Print_Real_Unit('theta_k_2',theta_k_2,'degrees')
Print_Real('dk_2_lambda_g0',dk_2_lambda_g0)
Print_Real_Unit('dk2',dk2,'m')

Print_Real_Unit('theta_k_3',theta_k_3,'degrees')
Print_Real('dk_3_lambda_g0',dk_3_lambda_g0)
Print_Real_Unit('dk3',dk3,'m')

Print_Real_Unit('theta_k_4',theta_k_4,'degrees')
Print_Real('dk_4_lambda_g0',dk_4_lambda_g0)
Print_Real_Unit('dk4',dk4,'m')

Print_Real_Unit('theta_k_5',theta_k_5,'degrees')
Print_Real('dk_5_lambda_g0',dk_5_lambda_g0)
Print_Real_Unit('dk5',dk5,'m')

Print_Real_Unit('z1',z1,'m')
Print_Real_Unit('z2',z2,'m')
Print_Real_Unit('z3',z3,'m')
Print_Real_Unit('z4',z4,'m')
Print_Real_Unit('z5',z5,'m')
Print_Real_Unit('z6',z6,'m')
Print_Real_Unit('z7',z7,'m')

Print_Real_Unit('L',L,'m')
Print_Real_Unit('dp',dp,'m') 

 IFigure = IFigure + 1;
 figure_max(IFigure)
 plot(a_g_HFSS_Plot/m, K_Z0_HFSS_Plot, 'r', 'linewidth', 9);
 hold on
 plot(a_g/m, K_Z0, 'bo', 'linewidth', 13);
 hold off
 title('Cubic Spline Interpolation')
 % xlabel('{\ita_g}  (mm)')
 % ylabel(' {\itK} / {\itZ}_{0}   (\Omega/\Omega)',...
 %    
 %'VerticalAlignment', 'bottom')
 % legend(...
 %     ' CS Interpolation',...
 %     ' WG Iris Design', ...
 %     'Location', 'Best')
 grid on
 grid minor
 set(gca, 'FontName', 'times new roman', 'FontSize', NF)
 set(gca, 'linewidth', 2.5)


 %--------------------------------------------------------------------
IFigure = IFigure + 1;
figure_max(IFigure)
plot(a_g_HFSS_Plot/m, phi_Plot, 'r', 'linewidth', 9);
hold on
plot(a_g/m, phi, 'bo', 'linewidth', 13);
hold off
title('Cubic Spline Interpolation')
% xlabel('{\ita_g} (mm)')
% ylabel(' {\it\phi} (°)', 'VerticalAlignment', 'bottom')
% legend(...
% ' CS Interpolation',...
% ' WG Iris Design', ...
% 'Location', 'Best')
grid on
grid minor
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
set(gca, 'linewidth', 2.5)