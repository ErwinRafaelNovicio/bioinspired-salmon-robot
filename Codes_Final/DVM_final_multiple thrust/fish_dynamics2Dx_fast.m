function [M,b,Torque_al] = fish_dynamics2Dx_fast(x,t,m_vec,l_vec,J_vec,h,C_am,C_d,A_vec,Theta_ref,kp,F_vort,M_vort,N_tail)
%% FISH_DYNAMICS2DX Bioinspired Robotics
% prof. S. Cinquemani, ing. G. Bianchi
% A.Y. 2025/2026
%% 21-10-2025 - 2D Discrete Vortex Model (DVM) simulation of a swimming BCF fish
% 
% 
% This function takes as input the state of the system, its geometry, and the 
% single forces appiled.
% 
% This function provides as output the mass matrix and the vector containing 
% all the Lagrangian components.
n_link = length(m_vec);
rho = 1000;                % Density of the fluid [kg/m^3]
%% 
% If the added mass and drag coefficient are scalar, it is assumed that they 
% are equal for all the links.
if isscalar(C_am)
    C_am = C_am*ones(1,n_link);
    C_am(end) = 0;
end
m_vec = (1+C_am).*m_vec;
if isscalar(C_d)
    C_d = C_d*ones(n_link,1);
end
%% 
% Extraction of position and velocities from the state
theta = x(3:n_link+2);
omega = x(n_link+5:2*n_link+4);
% Mass matrix
% 
M = zeros(n_link+2,n_link+2);
M(1:2,1:2) = sum(m_vec)*eye(2);
for ii = 1:n_link
    M(1,ii+2) = -l_vec(ii)*sin(theta(ii))*(sum(m_vec(ii:end))-m_vec(ii)/2);
    M(ii+2,1) = M(1,ii+2);
    M(2,ii+2) = l_vec(ii)*cos(theta(ii))*(sum(m_vec(ii:end))-m_vec(ii)/2);
    M(ii+2,2) = M(2,ii+2);
end
for ii = 1:n_link
    for jj = 1:n_link
        if ii == jj
            M(2+ii,2+jj) = J_vec(ii) + l_vec(ii)^2*(sum(m_vec(ii:end))-m_vec(ii)*3/4);
        else
            M(2+ii,2+jj) = l_vec(ii)*l_vec(jj)*cos(theta(ii)-theta(jj))*(sum(m_vec(max(ii,jj):end))-m_vec(max(ii,jj))/2);
        end
    end
end
% Centrifugal and Coriolis terms
% 
F_in_cent = zeros(n_link+2,1);
for ii  = 1:n_link
    F_in_cent(1) = F_in_cent(1) - omega(ii)^2*l_vec(ii)*cos(theta(ii))*(sum(m_vec(ii:end))-m_vec(ii)/2);
    F_in_cent(2) = F_in_cent(2) - omega(ii)^2*l_vec(ii)*sin(theta(ii))*(sum(m_vec(ii:end))-m_vec(ii)/2);
end
for ii = 1:n_link
    for jj = 1:n_link
        F_in_cent(2+jj) = F_in_cent(2+jj) + l_vec(ii)/2*l_vec(jj)*omega(jj)^2*sin(theta(ii)-theta(jj))*(sum(m_vec(jj:end))-m_vec(jj)/2);
    end
end
% Hydrodynamic drag forces
% $$F_{{\;}_{\textrm{drag}} } =-\frac{1}{2}\rho \;c_d \;A\;v_t \;|v_t |$$
vo0x = x(n_link+3);
vo0y = x(n_link+4);
vt =  vo0x*cos(theta(1))+vo0y*sin(theta(1));
vn = -vo0x*sin(theta(1))+vo0y*cos(theta(1));
F_drag_t = -0.5*rho*C_d(1)*A_vec(1)*vt*abs(vt);
F_drag_n = -0.5*rho*(2.4*l_vec(ii)*h)*vn*abs(vn);
F_drag_X = F_drag_t*cos(theta(1))-F_drag_n*sin(theta(1));
F_drag_Y = F_drag_t*sin(theta(1))+F_drag_n*cos(theta(1));
M_drag = -0.5*rho*C_d.*A_vec'.*omega.*abs(omega).*l_vec'/2;
F_drag = [F_drag_X; F_drag_Y; M_drag];
if any(abs(F_drag)>10^3)
    warning('t=%.3e: Unexpected drag force F_drag=%.2e N', t, max(abs(F_drag)));
end
% Applied torques
% $M=k_p \left(\theta_{\textrm{ref}} \left(t\right)\right)-\theta \;\left(t\right)$)
theta_rel = zeros(1,n_link);
Torque_c = zeros(1,n_link);
for ii = 2:n_link
    theta_rel(ii) = theta(ii)-theta(ii-1);
    Torque_c(ii) = kp(ii)*(Theta_ref(ii)-theta_rel(ii));
end
%% 
Torque_al=Torque_c;
% The torque is multiplied by a sigmoid function to avoid abrupt starting
Torque_c = Torque_c*(-1+2./(1+exp(-100*t)));
F_torque = zeros(n_link+2,1);
for ii = 1:n_link-1
    F_torque(2+ii) = Torque_c(ii)-Torque_c(ii+1);
end
F_torque(2+n_link) = Torque_c(n_link);
% Hydrodynamic vortex-induced forces
% 
F_vortex = zeros(n_link+2,1);
F_vortex(1) = sum(F_vort(1,:));
F_vortex(2) = sum(F_vort(2,:));
%N_tail=5;
N_rigid=n_link-N_tail;
for I_ele=1:N_tail
for ii = 1:n_link-1
    F_vortex(2+ii) = -F_vort(1,I_ele)*l_vec(ii)*sin(theta(ii))+F_vort(2,I_ele)*l_vec(ii)*cos(theta(ii));
end
F_vortex(2+N_rigid+I_ele) = M_vort(I_ele);
end
% Generalized force vector
b = -F_in_cent + F_drag + F_torque + F_vortex;
end