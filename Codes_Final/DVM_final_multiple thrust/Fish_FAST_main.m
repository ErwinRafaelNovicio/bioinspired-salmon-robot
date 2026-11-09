%% Bioinspired Robotics
% prof. S. Cinquemani, ing. G. Bianchi
% A.Y. 2025/2026
%% 21-10-2025 - 2D Discrete Vortex Model (DVM) simulation of a swimming BCF fish (SUBCARANGIFORM ADAPTED)
%% Main script

clear all
close all
clc

addpath("funzioniDVM\")

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

%% Select if you want to RUN the simulation or only look to the results of a previous simulation
% Put the flag on the 'RunSimulation' and 'SaveSimulation' variables and 
% write the title of the file that you want ot save or to load
RunSimulation = true;
SaveSimulation = false;% se vera salva vortices, il che può pesare qualche giga!
FileName = 'Sim_3s.mat';

% Decide also if you want to save a gif of the animation
save_gif = false; % Flag to decide if save or not the gif, put "true" or "false"
gifname = 'animation.gif';

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
l_vec = [10,    10,   30,   10,    10,    10,    4,    4,   4,    4,    4]/100 * l_tot; % Length of the links [m]
J_vec = 1/12 * m_vec .* l_vec.^2; % Moment of inertia of the links [kg*m^2]
h = 8.3 * 1e-3; % Thickness of the fin and body [m]

% A_vec = [1565, 3856, 15620, 5050, 4640, 3888, 1270, 1065, 843, 516, 172]*1e-6; % Cross-section area of the links [m^2]
A_vec = [1565, 3856, 15620, 5050, 4640, 3888, 1270, 1065, 843, 516, 2000]*1e-6; % Cross-section area of the links [m^2]

C_am = 0.1;                     % Added mass coefficient
C_d = 0.1;                      % Drag coefficient

% Stiffnesses
k_flex = 5e-3;
k_vec = 0.8*[ 1.9      1.55     0.44      0.2848    0.077006  0.077006       0.015      0.0075     0.00375     0.00165      0.0015]; % Torsional stiffness of the joints [Nm/rad]
figure
bar(cumsum(l_vec),k_vec)
hold on
plot([0,cumsum(l_vec)],zeros(n_link+1,1),'o',LineWidth=2)
plot([0,l_tot],[0,0],'k','LineWidth',2)

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


M = zeros(n_link+2,n_link+2);
M(1:2,1:2) = sum(m_vec)*eye(2);
for ii = 1:n_link
    M(1,ii+2) = -l_vec(ii)*sin(0)*(sum(m_vec(ii:end))-m_vec(ii)/2);
    M(ii+2,1) = M(1,ii+2);
    M(2,ii+2) = l_vec(ii)*cos(0)*(sum(m_vec(ii:end))-m_vec(ii)/2);
    M(ii+2,2) = M(2,ii+2);
end
for ii = 1:n_link
    for jj = 1:n_link
        if ii == jj
            M(2+ii,2+jj) = J_vec(ii) + l_vec(ii)^2*(sum(m_vec(ii:end))-m_vec(ii)*3/4);
        else
            M(2+ii,2+jj) = l_vec(ii)*l_vec(jj)*cos(0)*(sum(m_vec(max(ii,jj):end))-m_vec(max(ii,jj))/2);
        end
    end
end
[modes,lam]=eig(M\K);
lam=diag(lam)/2/pi; % mechanical eigenfreq

%% Define the actuator generalized force 




Amp_vec=[0,0,0,0.0347    0.0782    0.1515];%height=25 in cinematica
freq = 2  ; % Frequency [Hz] 1.75
omega = 2 * pi * freq;
phi = 20 / 180 * pi; % Phase shift




switch link_active
    case 'head' % all body is actuated

        ind = [1, 2, 4 : 6]; % Index of the actuated link

        Theta_ref_fun = @(t) [Amp_vec(1) * sin(omega*t); ...
                              Amp_vec(2) * sin(omega*t); ...
                              0;
                              Amp_vec(4) * sin(omega*t); ...
                              Amp_vec(5) * sin(omega*t); ...
                              Amp_vec(6) * sin(omega*t); ...
                              0; 0; 0; 0; 0];   % last two joints are passive
    case 'tail' %only the tail is actuated

        ind = 4 : 6; % Index of the actuated link

        Theta_ref_fun = @(t)                   [0; 0; 0;
                                               Amp_vec(4) * sin(omega*t); ...
                                               Amp_vec(5) * sin(omega*t+phi); ...
                                               Amp_vec(6) * sin(omega*t+2*phi); ...
                                               0; 0; 0; 0; 0];   % last five joints are passive
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

%% --> Simulation start
% If the Flag is true the simulation runs and save all the results in a
% matrix that can then be loaded to skip the passages if the parameters
% don't change

if RunSimulation == true

    % Initialization of the variables containing the result
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
    
    % Integrate numerically the equation of motion
    %keep track of the simulation progress
    tic;
    nSteps = length(tv);
    progress_marks = 0 : 10 : 100;
    
    Torque_active_links = zeros(nSteps,11);
    
    for ii = 2 : length(tv)
    
        vortices_old = vortices{ii-1};
        dxdt_old = dxdt;
        [dxdt, vortices{ii}, Torque_active_links(ii,:)] = fish_vortices_dyn_fast(tv(ii), x, h, m_vec, l_vec, J_vec, K, C_am, C_d, A_vec, Theta_ref_fun, kp, Nb, sigma, vortices_old, dt, T_fade,N_tail);
    
        % Trapezoidal integration rule
        x = x + (0.5 * dxdt_old + 0.5 * dxdt) * dt;
       
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

    if SaveSimulation == true

        % Saving the results in a matrix
        disp('')
        disp('')
        disp(['--> Saving the results in a file: ', FileName])
    
        save(FileName, 'vortices', 'X_hist', 'Torque_active_links', '-v7.3');
    
        disp('-----> Saving complete!')

    end

else % If the 'RunSimulation' is FALSE you have to load the file that contains the results of the simulation

    disp('')
    disp('')
    disp(['--> Loading the results from the file: ', FileName])

    load(FileName);

    disp('----> Loading done!')

end

%--------------------------------------------------------------------------
%% PLot of the results and saving a gif

delayTime = 0.05;   % secondi tra i frame
writeEvery = 50;    % salva ogni quanto (qui coincide con il passo del tuo for)
firstFrame = true;
follow= @(I) (tv(I)-min(tv(I),3)) * X_hist(14,ii); % sposta la "camera" dopo 3 secondi dell'animation
anime = figure;
set(gca, 'color', [100,200,255]/255)
tt = 0;
x_link = zeros(1,n_link);
y_link = zeros(1,n_link);

for ii = 1 : 50 : length(vortices)-1

    hold on; 
    axis equal;
    if ii > 1
        delete(p)
        delete(p_hinge)
        delete(p_v)
    end

    % --- Reconstruct link positions (come prima) ---
    theta_link = X_hist(3 : n_link+2 , ii);
    x0 = X_hist(1,ii);
    y0 = X_hist(2,ii);
    x_link(1) = x0 + l_vec(1)*cos(theta_link(1));
    y_link(1) = y0 + l_vec(1)*sin(theta_link(1));

    for jj = 2 : n_link
        x_link(jj) = x_link(jj-1) + l_vec(jj)*cos(theta_link(jj));
        y_link(jj) = y_link(jj-1) + l_vec(jj-1)*sin(theta_link(jj));
    end

    % --- Draw links and hinges ---
    p(1) = plot([x0 x_link(1)], [y0 y_link(1)], 'k-', 'LineWidth', 4);

    for jj = 2:n_link
        p(jj) = plot([x_link(jj-1) x_link(jj)], [y_link(jj-1) y_link(jj)], 'k-', 'LineWidth', 3);
    end

    p_hinge = plot([x0 x_link], [y0 y_link], 'ro', 'MarkerFaceColor', 'r');

    % --- Draw vortices with a custom color map ---
    gk = vortices{ii}.gamma;
    pk = vortices{ii}.pos;
    p_v = scatter(pk(2:end,1), pk(2:end,2), 20, gk(2:end), 'filled');
    colormap(blue_white_red(256)); 
    colorbar;
    clim([-5e-5, 5e-5])
    title(sprintf('t = %.2f s', tv(ii)));
    xlim([-0.3 1.5*sum(l_vec)]+follow(ii)); ylim([-0.1 0.1]); 
    xlabel('x'); ylabel('y');
    drawnow

    % --- Capture and write GIF ---
    if mod(ii-1, writeEvery) == 0 && save_gif
        frame = getframe(anime);                 % cattura figura
        img = frame2im(frame);                   % rgb image
        [imind, cm] = rgb2ind(img, 256);         % indexed image per GIF
        if firstFrame
            imwrite(imind, cm, gifname, 'gif', 'LoopCount', inf, 'DelayTime', delayTime);
            firstFrame = false;
        else
            imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
        end
    end

end

%% Plot relevant data as functon of time

color_ln = ['#0072BD'; '#D95319'; '#EDB120'; '#7E2F8E'; '#77AC30'; '#4DBEEA'; '#A2142F'; '#008080'; '#B300B3'; '#996600'; '#800000'];

n_active = length(ind);
theta_absolute = X_hist(3 : n_link+2, :); % Joint angles [rad]
theta_relative = [theta_absolute(1, :); diff(theta_absolute, 1, 1)]; % θ_rel(i) = θ_(i+1) - θ_i

theta_tail = theta_relative(end, :); % Relative angle of the last link of the felxible tail
omega_tail = X_hist(2*(2+n_link),:) - X_hist(2*(2+n_link)-1,:);

% Plot of the displacement in the x, y direction as function of time and of
% the relative angle of the last link (CAUDAL FIN)
figure;
plot(tv, X_hist(1,:)) % Displacement in the X direction
hold on
plot(tv, X_hist(2,:)) % Displacement in the Y direction
hold on
plot(tv, theta_tail)
xlabel('Time [s]')
title('Displacements and relative angle tail')
legend('x', 'y', '$\theta_{rel}$')

figure;
plot(X_hist(1, :), X_hist(2, :))
xlabel('X displacement [m]')
ylabel('Y displacement [m]')
title('trajectory in the XY plane')
  
% Plot of the velocities in the X and Y direction of the all fish as
% function of time and of the relative angular velocity of the last link
figure;
plot(tv, X_hist(1+n_link+2,:)) % Velocity in the X direction
hold on
plot(tv, X_hist(2+n_link+2,:)) % Velocity in the Y direction
hold on
plot(tv, omega_tail)
xlabel('Time [s]')
title('Velocities and relative anglular velocity tail')
legend('$v_x$', '$v_y$', '$\omega_{rel}$')

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

% Plot of the normal and tangential force for each passive link
figure;
ax1 = axes;
plot(ax1, tv(1:length(vortices)), F_vort_t, 'LineStyle', '-')
xlabel('Time [s]')
ylabel('Hydrodynamic force [N]')
title('Normal and tangential forces')
hold(ax1, 'on')
h1 = plot(ax1, tv(1:length(vortices)), F_vort_n, 'LineStyle', '-.', 'HandleVisibility', 'off');
for cont = 1 : N_tail
    h1(cont).Color = color_ln(cont, :);
end

legend(ax1, '$F_{t-7}$', '$F_{t-8}$', '$F_{t-9}$', '$F_{t-10}$', '$F_{t-11}$', 'Location', 'northeastoutside')

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % For the second legend (NEGLECT IT)
    ax2 = axes('Color', 'none', 'XAxisLocation', 'bottom', 'YAxisLocation', 'left', 'Visible', 'off');
    hold(ax2, 'on')
    plot(ax2, NaN, NaN * ones(N_tail, 1), 'LineStyle', '-.')
    legend(ax2, '$F_{n-7}$', '$F_{n-8}$', '$F_{n-9}$', '$F_{n-10}$', '$F_{n-11}$', 'Location', 'southeastoutside', 'Color', 'white')
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Plot of the forces in the X and Y direction for each passive link
figure;
ax1 = axes;
plot(ax1, tv(1:length(vortices)), F_vort_x, 'LineStyle', '-')
xlabel('Time [s]')
ylabel('Hydrodynamic force [N]')
title('X and Y direction forces')
hold(ax1, 'on')
h1 = plot(ax1, tv(1:length(vortices)), F_vort_y, 'LineStyle', '-.');
for cont = 1 : N_tail
    h1(cont).Color = color_ln(cont, :);
end

legend(ax1, '$F_{x-7}$', '$F_{x-8}$', '$F_{x-9}$', '$F_{x-10}$', '$F_{x-11}$', 'Location', 'northeastoutside')

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % For the second legend (NEGLECT IT)
    ax2 = axes('Color', 'none', 'XAxisLocation', 'bottom', 'YAxisLocation', 'left', 'Visible', 'off');
    hold(ax2, 'on')
    plot(ax2, NaN, NaN * ones(N_tail, 1), 'LineStyle', '-.')
    legend(ax2, '$F_{y-7}$', '$F_{y-8}$', '$F_{y-9}$', '$F_{y-10}$', '$F_{y-11}$', 'Location', 'southeastoutside', 'Color', 'white')
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

Gamma = vortices{end}.gammaTotal;

figure
plot(tv(1:length(vortices)),Gamma)
xlabel('Time [s]')
ylabel('$\Gamma$ [m$^2$/s]')
title('Total circulation $\Gamma$')
ylim([-1e-5 1e-5])

%% Extract time and joint angles
time = tv;

% Plot joint angle evolution
figure;
subplot(1, 2, 1)
h1 = plot(time, theta_absolute(ind, :)'); % for the active joints
grid on
xlabel('Time [s]');
ylabel('Joint angles [rad]');
title('Joint absolute angles of active links');

subplot(1, 2, 2)
h2 = plot(time, theta_absolute(~ismember(1:end, ind), :)');
grid on
xlabel('Time [s]');
ylabel('Joint angles [rad]');
title('Joint absolute angles of the passive links');
for cont = n_active+1 : n_link
    h2(cont-n_active).Color = color_ln(cont, :);
end

switch link_active
    case 'head'
        
        legend(h1, '$\theta_1$', '$\theta_2$', '$\theta_4$', '$\theta_5$', '$\theta_6$', 'Interpreter','latex')
        legend(h2, '$\theta_3$', '$\theta_7$', '$\theta_8$', '$\theta_9$', '$\theta_{10}$', '$\theta_{11}$', 'Interpreter','latex')

    case 'tail'

        legend(h1, '$\theta_4$', '$\theta_5$', '$\theta_6$', 'Interpreter','latex')
        legend(h2, '$\theta_1$', '$\theta_2$', '$\theta_3$', '$\theta_7$', '$\theta_8$', '$\theta_9$', '$\theta_{10}$', '$\theta_{11}$', 'Interpreter','latex')

end

% Plot joint angle evolution (relative)
figure;
subplot(1, 2, 1)
h1 = plot(time, theta_relative(ind, :)'); % for the active joints
grid on
xlabel('Time [s]');
ylabel('Joint angles [rad]');
title('Joint relative angles of active links');

subplot(1, 2, 2)
h2 = plot(time, theta_relative(~ismember(1:end, ind), :)');
grid on
xlabel('Time [s]');
ylabel('Joint angles [rad]');
title('Joint relative angles of the passive links');
for cont = n_active+1 : n_link
    h2(cont-n_active).Color = color_ln(cont, :);
end

switch link_active
    case 'head'
        
        legend(h1, '$\theta_1$', '$\theta_2$', '$\theta_4$', '$\theta_5$', '$\theta_6$', 'Interpreter','latex')
        legend(h2, '$\theta_3$', '$\theta_7$', '$\theta_8$', '$\theta_9$', '$\theta_{10}$', '$\theta_{11}$', 'Interpreter','latex')

    case 'tail'

        legend(h1, '$\theta_4$', '$\theta_5$', '$\theta_6$', 'Interpreter','latex')
        legend(h2, '$\theta_1$', '$\theta_2$', '$\theta_3$', '$\theta_7$', '$\theta_8$', '$\theta_9$', '$\theta_{10}$', '$\theta_{11}$', 'Interpreter','latex')

end

%% Amplitude of the active part of the tail of the fish
x_active = zeros(1, n_active + 1);
y_active = zeros(1, n_active + 1);

% Find the value at which the fish has the maximum amplitude
[~, idx] = max(theta_absolute(4, :));

figure()

for cont = 1 : n_active

    theta_max(cont) = theta_absolute(cont + n_active, idx);

    x_active(cont + 1) = x_active(cont) + l_vec(cont + n_active) * cos(theta_max(cont));
    y_active(cont + 1) = y_active(cont) + l_vec(cont + n_active) * sin(theta_max(cont));

    textString = sprintf('$\\theta_{abs} = %.2f^{\\circ}$', rad2deg(theta_max(cont)));
    text(x_active(cont + 1), y_active(cont + 1) + 0.0005, textString, 'Interpreter', 'latex', 'FontSize', 10, 'HorizontalAlignment', 'right')
    hold on

end

plot(x_active, y_active, 'LineWidth', 2)
grid on
xlabel('x [m]', 'Interpreter', 'latex')
ylabel('y [m]', 'Interpreter', 'latex')
title('Amplitude of the active links of the tail', 'Interpreter', 'latex', 'FontWeight', 'bold')
hold on
scatter(x_active, y_active, 'filled', 'r')

%% coppia motrice
figure

plot(tv,Torque_active_links(:,ind))
title('Torque @active link due to motor')
legend('link 4','link 5','link 6')
%save('motore.mat',"tv","motorTorque")
T_fish=rms(Torque_active_links((tv(1:end-1)>3/freq),ind)); %da inserire in cinematica

%% Thrust
% Compute the instantaneous thrust for each time instant summing the force
% in the X direction along the rows
thrust_inst = sum(F_vort_x, 2);

% Identify the Steady-State Window exluding the transient part, so starting
% from the second period
t_start = 0.5; 
t_end = tv(end);
steady_idx = find(tv >= t_start & tv <= t_end);

% Numerical Integration for Mean Thrust
mean_thrust = trapz(tv(steady_idx), thrust_inst(steady_idx)) / (tv(steady_idx(end)) - tv(steady_idx(1)));

% Display Results
fprintf('Mean Thrust Generated: %.5f N\n', mean_thrust);

% Plotting results
figure()
plot(tv, thrust_inst, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Instantaneous')
hold on
plot(tv(steady_idx), thrust_inst(steady_idx), 'b', 'LineWidth', 1.5, 'DisplayName', 'Steady-State')
yline(mean_thrust, 'r', 'LineWidth', 2, 'DisplayName', 'Mean Thrust');
xlabel('Time [s]');
ylabel('Force [N]');
title('Total Thrust');
legend;
grid on;