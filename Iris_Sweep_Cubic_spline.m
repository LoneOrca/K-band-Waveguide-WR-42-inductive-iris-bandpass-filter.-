clc;
clear all;
close all;
clf;
%----------------------------------------------------------------------
%% Iris Sweep Values - Cubic spline Interpolation 
%----------------------------------------------------------------------
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

 
Xs_Z0 = ppval(P_Xs_Z0_vs_ag, a_g); %%%
 Xp_Z0 = ppval(P_Xp_Z0_vs_ag, a_g); %%%
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
 
