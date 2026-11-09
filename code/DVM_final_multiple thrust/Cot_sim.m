%% Bioinspired Robotics
% prof. S. Cinquemani, ing. G. Bianchi
% A.Y. 2025/2026
%% 21-10-2025 - 2D Discrete Vortex Model (DVM) simulation of a swimming BCF fish (SUBCARANGIFORM ADAPTED)
%% Main script
clear all
close all
clc

amp_mat=[
0.008623682469301	0.019619734666751	0.039124498466539;
0.021360918519338	0.048401733417289	0.095088573706081;
0.026609208804730	0.060196452330995	0.117583781659106;
0.042962384279214	0.096720996515386	0.185834210015383;
0.048641894672386	0.109329845789560	0.208948885354694;
0.054450801966385	0.122187389779294	0.232306428074812];%in cinematica corrispondo a height:[10,15,20,22,28,30,32];

frequencies=[1,2,5,7,10,12];

heights=[15,20,22,28,30,32];

Tq=[
2.445190535319982	1.386019584470329	0.201072951057538;
2.444193519411294	1.387027663112475	0.212405961848426;
2.443529861442145	1.387725458374725	0.219855967618202;
2.440492708391480	1.391211710968531	0.253296109699718;
2.439088685110279	1.392999363933096	0.268534044984165;
2.437463575773767	1.395222158851176	0.286103218111153];

%data storage
dist=zeros([size(amp_mat,1),length(frequencies)]);
T_fish=zeros([size(amp_mat,1),length(frequencies),3]);
Cm_lag=zeros([size(amp_mat,1),length(frequencies)]);
Pm_lag=zeros([size(amp_mat,1),length(frequencies)]);
CoT=zeros([size(amp_mat,1),length(frequencies)]);

for Jcota=1:size(amp_mat,1)
for Jcotf=1:length(frequencies)
 clearvars('-except','Jcotf','Jcota','T_fish','frequencies','amp_mat','Tq','CoT','dist','Cm_lag','Pm_lag','heights')  


set(0,'DefaultAxesFontSize',18);
set(0,'DefaultTextFontSize',18);
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 0.8);
set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1);
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesXGrid','on');
set(0,'DefaultAxesYGrid','on');
set(0,'DefaultAxesGridLineStyle','--');
set(0,'DefaultLineMarkerSize',8);
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Define the fish as a chain of rigid links connected by hinges
% The fish is swimming leftwards, in the negative X direction

n_link = 11;
N_tail = n_link - 6; % number of element needed to discretize flexible tail, note that n_rigid_el = 6

% Select if the two links in the head are moved or not
link_active = 'tail'; % put what links you want to activate
                      % 'head' also the two links in the head are activated
                      % 'tail' only the links in the tail are moving

% Computation of the mass vector
rho_fish = 997; % [kg/m3]
V_fish = 0.00154451; % [m3]
m_tot = rho_fish * V_fish; % [kg]
m_vec = [2.3, 5.7, 70, 6, 6, 4.4, 1.3, 1.2, 1.1, 1,1]/100 * (m_tot); % Mass of the links [kg]

% Definition of the geometry of the fish
l_tot = 0.5; % [m]
l_vec = [10,    10,   30,   10,    10,    10,    4,    4,   4,    4,    4]/100*l_tot; % Length of the links [m]
J_vec = 1/12 * m_vec .* l_vec.^2; % Moment of inertia of the links [kg*m^2]

h = 8.3 * 1e-3; % Thickness of the fin and body [m]
A_vec = [1565, 3856, 15620, 5050, 4640, 3888, 1270, 1065, 843, 516, 2500]*1e-6; % Cross-section area of the links [m^2]

C_am = 0.1;                     % Added mass coefficient
C_d = 0.1;                      % Drag coefficient
k_flex = 15e-4;

k_vec = 0.8*[ 1.9      1.55     0.44      0.2848    0.077006  0.077006       0.015      0.0075     0.00375     0.00165      0.0015]; % Torsional stiffness of the joints [Nm/rad]




%% Write the equation of motion using Lagrange equation  
% Precompute stiffness matrix

K_mat = zeros(n_link,n_link);
K_mat(1,1) = k_vec(2);
K_mat(1,2) = -k_vec(2);
for ii = 2:n_link-1
    K_mat(ii,ii-1) = -k_vec(ii);
    K_mat(ii,ii) = k_vec(ii)+k_vec(ii+1);
    K_mat(ii,ii+1) = -k_vec(ii+1);
end
K_mat(n_link,n_link-1) = -k_vec(end);
K_mat(n_link,n_link) = k_vec(end);  

K = [zeros(2,2),       zeros(2,n_link);
     zeros(n_link,2),           K_mat];



%% Define the actuator generalized force 


Amp_vec=[0,0,0,amp_mat(Jcota,:)];
freq =frequencies(Jcotf)  ; % Frequency [Hz] 1.75
omega = 2 * pi * freq;
phi=20/180*pi;

jota=0;
transient=@(t,a) 1-exp(log(0.01)/(a+eps)*(t+a/100));%a è l'istante di tempo dove transient(t=a,a)=0.99
switch link_active
    case 'head'

        ind = [1, 2, 4 : 6]; % Index of the actuated link

        Theta_ref_fun = @(t) [Amp_vec(1) * sin(omega*t); ...
                              Amp_vec(2) * sin(omega*t); ...
                              0;
                              Amp_vec(4) * sin(omega*t); ...
                              Amp_vec(5) * sin(omega*t); ...
                              Amp_vec(6) * sin(omega*t); ...
                              0; 0; 0; 0; 0];   % last two joints are passive
    case 'tail'

        ind = 4 : 6; % Index of the actuated link

        Theta_ref_fun = @(t) transient(t,0)*[0; 0; 0;
                              Amp_vec(4) * sin(omega*t+jota); ...
                              Amp_vec(5) * sin(omega*t+phi+jota); ...
                              Amp_vec(6) * sin(omega*t+2*phi+jota); ...
                              0; 0; 0; 0; 0];   % last two joints are passive
end


% Theta_ref_fun = @(t) angoli_robot(t,omega,n_link,ind,35);
% ind = 4 : 6; % Index of the actuated link
kp = zeros(1, n_link);
kp(ind) = 1; % Proportional gain of the control

% Define the DVM parameters
Nb = 10; % Number of panels  per  link

%% Regularization function: $g_{\sigma \;} \left(x\right)=1-e^{-\frac{x^2 }{2\sigma{\;}^2 }}$

sigma = 0.1*l_vec(end)/Nb;    % Radius of the regularized vortex core [m]

%% Fading time after which the vortex dissolves, simulating the effect of viscosity.
% For T_fade != infty the total circulation is not conserved.

T_fade = Inf;

% Rearrange the equation of motion
% Define the initial conditions
x_init = zeros(2*(n_link+2),1);
x_init(n_link+3) = -0.25; % Initial velocity [m/s]

% -------------------------------------------------------------------------
Tmax = 3; % Simulation time [s]
% -------------------------------------------------------------------------

dt = 8e-4;
tv = 0 : dt : Tmax;

%% Initialization of the variables containing the result
X_hist = zeros(length(x_init), length(tv));
X_hist(:,1) = x_init;

vortices = cell(1,length(tv));
vortices{1}.pos = [0,0];
vortices{1}.gamma = 0;
vortices{1}.gamma_b = [];
vortices{1}.age = 0;
vortices{1}.F_vort = [0, 0];
vortices{1}.gammaTotal = 0;
vortices{1}.dGammadt = [];

dxdt = zeros(size(x_init));
x = x_init;





%% Integrate numerically the equation of motion

 if ~(freq>10 && heights(Jcota)>29.9)%non simulo con dati che mi portano a fallire la dvm
tic;

nSteps = length(tv);
progress_marks = 0 : 10 : 100;
Torque_active_links=zeros(nSteps,11);

for ii = 2 : length(tv)

    vortices_old = vortices{ii-1};
    dxdt_old = dxdt;
   [dxdt, vortices{ii},Torque_active_links(ii,:)] = fish_vortices_dyn_fast(tv(ii), x, h, m_vec, l_vec, J_vec, K, C_am, C_d, A_vec, Theta_ref_fun, kp, Nb, sigma, vortices_old, dt, T_fade,N_tail);

    % Trapezoidal integration rule
    x = x + (0.3 * dxdt_old + 0.7 * dxdt) * dt;
   
    X_hist(:,ii) = x; 

    % X_hist contains in the order:
    % - 1:2 = X and Y displacement
    % - 3:13 = absolute angles of the links of the fish
    % - 14:15 = X and Y velocities of the fish
    % - 16:26 = Angular velocities of the links of the fish

    progress = 100 * ii / nSteps;

    if any(abs(progress - progress_marks) < 0.5/nSteps*100)  % print ~every 10%
        elapsed = toc;
        disp(['Progress: ', num2str(progress, '%.0f'), '% completed (', ...
              num2str(elapsed, '%.1f'), ' s elapsed)']);
    end

end

disp('end dynamics')

% -------------------------------------------------------------------------

%--------------------------------------------------------------------------


%% Plot relevant data as functon of time

color_ln = ['#0072BD'; '#D95319'; '#EDB120'; '#7E2F8E'; '#77AC30'; '#4DBEEA'; '#A2142F'; '#008080'; '#B300B3'; '#996600'; '#800000'];

n_active = length(ind);
theta_absolute = X_hist(3 : n_link+2, :); % Joint angles [rad]
theta_relative = [theta_absolute(1, :); diff(theta_absolute, 1, 1)]; % θ_rel(i) = θ_(i+1) - θ_i

theta_tail = theta_relative(end, :); % Relative angle of the last link of the felxible tail
omega_tail = X_hist(2*(2+n_link),:) - X_hist(2*(2+n_link)-1,:);


% Computation and plot of the forces for each link that generates the thrust
F_vort_x = zeros(length(vortices), N_tail);
F_vort_y = zeros(length(vortices), N_tail);
F_vort_t = zeros(length(vortices), N_tail);
F_vort_n = zeros(length(vortices), N_tail);

% For each column there is the vortex force in direction X or Y applied on
% each flexible link
for cont = 1 : N_tail
    for ii = 1 : length(vortices)
        if ii == 1
            F_vort_x(ii, cont) = vortices{ii}.F_vort(end, 1);
            F_vort_y(ii, cont) = vortices{ii}.F_vort(end, 2);
        else 
            F_vort_x(ii, cont) = vortices{ii}.F_vort(end-N_tail+cont, 1);
            F_vort_y(ii, cont) = vortices{ii}.F_vort(end-N_tail+cont, 2);
        end
    end
    
    % The tangential and normal force have to be computed for each passive
    % link with respect its relative angle
    F_vort_t(:, cont) =  F_vort_x(:, cont) .* cos(theta_relative(cont, 1:length(vortices)))' + F_vort_y(:, cont) .* sin(theta_relative(cont, 1:length(vortices)))';
    F_vort_n(:, cont) = -F_vort_x(:, cont) .* sin(theta_relative(cont, 1:length(vortices)))' + F_vort_y(:, cont) .* cos(theta_relative(cont, 1:length(vortices)))';

end


Gamma = vortices{end}.gammaTotal;


%% Extract time and joint angles
time = tv;



%% efficienza
dist(Jcota,Jcotf)=X_hist(1,end)
T_fish(Jcota,Jcotf,:)=rms(Torque_active_links((tv(1:end-1)>min(3/freq,1.5)),ind))

Cm_lag(Jcota,Jcotf)=mean(squeeze(T_fish(Jcota,Jcotf,:))'./Tq(Jcota,:));
Pm_lag(Jcota,Jcotf)=Cm_lag(Jcota,Jcotf)*2*pi*omega;
CoT(Jcota,Jcotf)=Pm_lag(Jcota,Jcotf)*Tmax/(m_tot*X_hist(1,end))
else
dist(Jcota,Jcotf)=NaN
T_fish(Jcota,Jcotf,:)=[NaN,NaN,NaN];
CoT(Jcota,Jcotf)=NaN %plot con differenti ampiezze e freq

Cm_lag(Jcota,Jcotf)=NaN;
Pm_lag(Jcota,Jcotf)=NaN;
CoT(Jcota,Jcotf)=NaN
end

end
end
%%


save('data_CoT.mat','dist','T_fish','Cm_lag','Pm_lag','CoT','frequencies','heights')


%%
figure
surf(heights,frequencies,-CoT')
xlabel('amplitude')
ylabel('freq')
title('CoT')
figure
waterfall(heights,frequencies,-CoT')
xlabel('amplitude')
ylabel('freq')
title('CoT')