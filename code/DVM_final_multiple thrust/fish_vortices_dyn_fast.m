function [dxdt, vortices,Torque_c] = fish_vortices_dyn_fast(t, x, h, m_vec, l_vec, J_vec, K, C_am, C_d, A_vec, Theta_ref_fun, kp, Nb, sigma, vort_old, dt, T_fade,N_tail,Ti)
%% FISH_VORTICES_DYN Bioinspired Robotics
% prof. S. Cinquemani, ing. G. Bianchi
% A.Y. 2025/2026
%% 21-10-2025 - 2D Discrete Vortex Model (DVM) simulation of a swimming BCF fish
% 
% 
% This function takes as input the state of the system, computes the forces 
% and finds the acceleration.
% 
% The outputs of this function are the derivative of the state and a struct 
% containing the information about vortices at each time step.
n_link = length(m_vec);
Theta_ref = Theta_ref_fun(t);
%% 
% $$\mathit{\mathbf{x}}=\left\lbrack \begin{array}{c}x_0 \\y_0 \\\theta_1 \\\vdots 
% \\\theta_n \\\dot{x_0 } \\\dot{y_0 } \\\dot{\theta_0 } \\\vdots \\\dot{\theta_n 
% } \end{array}\right\rbrack \;\;\;\;\;\;\;\;\;\;$$
% Extraction of the positions and the velocities from the state vector
vO0 = x(n_link+3:n_link+4);           
pos = zeros(2, n_link+1);
vel = zeros(2, n_link+1);
vel(:,1) = vO0;
pos(:,1) = x(1:2);
omega = x(n_link+5:2*n_link+4);
theta = x(3:n_link+2);
%% 
% $$\left\lbrace \begin{array}{ll}x_i =x_{i-1} +l_{i\;} \cos \;\theta_{i\;}  
% & \;\\y_i =y_{i-1} +l_{i\;} \sin \;\theta_{i\;}  & \;\end{array}\right.\;\;\;\left\lbrace 
% \begin{array}{ll}\dot{x_i } ={\dot{x} }_{i-1} -{\dot{\theta} }_{i\;} l_{i\;} 
% \sin \;\theta_{i\;}  & \;\\\dot{y_i } ={\dot{y} }_{i-1} +{\dot{\theta} }_{i\;} 
% l_{i\;} \cos {\;\theta }_{i\;}  & \;\end{array}\right.$$
for ii = 1:n_link
    pos(:,ii+1) = pos(:,ii) + l_vec(ii)*[cos(theta(ii)); sin(theta(ii))];
    vel(:,ii+1) = vel(:,ii) + omega(ii)*l_vec(ii)*[-sin(theta(ii)); cos(theta(ii))];
end
N_rigid=n_link-N_tail;
xj = pos(1,N_rigid+1:end);
yj = pos(2,N_rigid+1:end);
vxj = vel(1,N_rigid+1:end);
vyj = vel(2,N_rigid+1:end);
thetaj = x((N_rigid+3):n_link+2);
omegaj = omega(N_rigid+1:end);
% Compute the hydrodynamic forces solving the vortex dynamics
l = l_vec(N_rigid+1:end);
[F_vort, M_vort, vort] = DVM_fast(t,l,xj,yj,thetaj,vxj,vyj,omegaj,Nb,sigma,h,vort_old,dt,T_fade,N_tail);
% Compute the mass matrix and the vector containing the inertial and the applied forces
[M, b,Torque_c] = fish_dynamics2Dx_fast(x,t,m_vec,l_vec,J_vec,h,C_am,C_d,A_vec,Theta_ref,kp,F_vort,M_vort,N_tail);
%% 
% Subtract the elastic forces from the right hand side
b = b-K*x(1:n_link+2);
% Compute the acceleration and the derivative of the state x
% Acceleration
a = M\b;         
if any(abs(a)>10^5)
    warning('t=%.3f: Unexpected acceleration a=%.2e m/s^2 or rad/s^2', t, max(abs(a)));
end
%% 
% State derivative
dxdt = zeros(length(x),1);
dxdt(1:n_link+2) = x(n_link+3:end);  % pos ← vel
dxdt(n_link+3:end) = a;              % vel ← acc
vortices = vort;
end