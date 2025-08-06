clc;
clear all;
close all;
clf;
%------------------------------------------------------------
%% EE-414-514 sn620094 Design Project 3
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
%fH = 22.5062*G;
%fL = 21.5062*G;



delta = ((fH - fL)/f0);

Print_Real_Unit('f0',f0,'Hz')
Print_Real_Unit('Bwf',BWf,'Hz')
Print_Real_Unit('fL',fL,'Hz')
Print_Real_Unit('fH',fH,'Hz')
Print_Real('delta',delta)
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
%
%----------------------------------------------------------------------
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
% %----------------------------------------------------------------------
% % Plots
% %----------------------------------------------------------------------
IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S11_Filter_dB, 'g', 'linewidth', 6)
hold on
plot(freq/G, S21_Filter_dB, 'r', 'linewidth', 6);
plot(f0/G,interp1(freq/G,S21_Filter_dB,f0/G),'ro','linewidth',10);
plot(fL/G,interp1(freq/G,S21_Filter_dB,fL/G),'ro','linewidth',10);
plot(fH/G,interp1(freq/G,S21_Filter_dB,fH/G),'ro','linewidth',10);
plot(fL/G,interp1(freq/G,S11_Filter_dB,fL/G),'go','linewidth',10);
plot(fH/G,interp1(freq/G,S11_Filter_dB,fH/G),'go','linewidth',10);

hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{k1} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
legend('|{\itS}_{11}|','|{\itS}_{21}|')
axis([f_min/G, f_max/G, -40, 0])
set(gca, 'xtick', f_min/G : 0.1 : f_max/G);
set(gca, 'linewidth', 2.0)


IFigure = IFigure + 1;
figure_max(IFigure)
plot(freq/G, S21_Filter_dB, 'r', 'linewidth', 6)
hold on
plot(f0/G,interp1(freq/G,S21_Filter_dB,f0/G),'ro','linewidth',10);
plot(fL/G,interp1(freq/G,S21_Filter_dB,fL/G),'ro','linewidth',10);
plot(fH/G,interp1(freq/G,S21_Filter_dB,fH/G),'ro','linewidth',10);
hold off
grid on
grid minor 
xlabel('{\itf} (GHz) ')
ylabel('| {\itS}_{21} | ( dB ) ', ...
'VerticalAlignment', 'bottom')
set(gca, 'FontName', 'times new roman', 'FontSize', NF)
legend('|{\itS}_{21}|')
axis([21.4*G/G, 22.6*G/G, -0.02, 0])
set(gca, 'xtick', f_min/G : 0.1 : f_max/G);
set(gca, 'linewidth', 2.0)
