function [F_vort, M_vort, vort] = DVM_fast(t,l,xj,yj,thetaj,vxj,vyj,omegaj,Nb,sigma,h,vort_old,dt,T_fade,N_tail)
%% DVM Bioinspired Robotics
% prof. S. Cinquemani, ing. G. Bianchi
% A.Y. 2025/2026
%% 21-10-2025 - 2D Discrete Vortex Model (DVM) simulation of a swimming BCF fish
% 
% 
% This function takes as input the statem the geometry and the DVM parameters.
% 
% This function provides as output the hydrodynamic force and moment and the 
% vortex positions and circulation.
rho = 1000;
%% 
% The previous vortices positions and circulations are stored
vortices = vort_old;
Gamma_b_old = vort_old.gamma_b;
if ~isempty(vortices.dGammadt)
    dGamma_dt_old = vortices.dGammadt(:,end)';
else
    dGamma_dt_old = zeros(Nb,1);
end
% if t <= dt && thetaj == 0
%     thetaj = 0.01;
% end
%% 
% The position and velocity of the link are stored
x_tail = xj(end) ;
y_tail = yj(end);
% Velocity in the collocation points (CP) 
% 
%N_tail=5;

fast_col='si';
switch fast_col
   case 'si'
Ex=diff(xj)/Nb.*ones(Nb,N_tail);
x=[xj(1:N_tail);xj(1:N_tail)+cumsum(Ex)];
X_cp=1/2*(x(1:end-1,:)+x(2:end,:));
x_cp=reshape(X_cp,1,[]);

Ey=diff(yj)/Nb.*ones(Nb,N_tail);
y=[yj(1:N_tail);yj(1:N_tail)+cumsum(Ey)];
Y_cp=1/2*(y(1:end-1,:)+y(2:end,:));
y_cp=reshape(Y_cp,1,[]);
dxp=reshape(Ex,1,[]);
dyp=reshape(Ey,1,[]);

Vcpx=vxj(1:N_tail)-(Y_cp-yj(1:N_tail)).*omegaj(:)';
Vcpy=vyj(1:N_tail)+(X_cp-xj(1:N_tail)).*omegaj(:)';
 v_cp_x=reshape(Vcpx,1,N_tail*Nb);
 v_cp_y=reshape(Vcpy,1,N_tail*Nb);
    case 'no'
dxp=[];
dyp=[];
x_cp=[];
y_cp=[];
v_cp_x=[];
v_cp_y=[];
for I_ele=1:N_tail

x_seg = linspace(xj(I_ele), xj(I_ele+1), Nb+1);
y_seg = linspace(yj(I_ele), yj(I_ele+1), Nb+1);

%% 
% CP are placed in the middle of each segment
x_cp =[x_cp, 0.5 * (x_seg(1:end-1) + x_seg(2:end))];
y_cp = [y_cp,0.5 * (y_seg(1:end-1) + y_seg(2:end))];
dxp =[dxp,(diff(x_seg))] ;
%disp(mean(diff(x_seg))-E(1,I_ele))
dyp = [dyp,diff(y_seg)];
%% 
% Velocity of the collocation points on the panel
v_cp_x = [v_cp_x,vxj(I_ele) - (0.5 * (y_seg(1:end-1) + y_seg(2:end)) - yj(I_ele)) * omegaj(I_ele)];
v_cp_y = [v_cp_y,vyj(I_ele) + (0.5 * (x_seg(1:end-1) + x_seg(2:end)) - xj(I_ele)) * omegaj(I_ele)];
end
% per verificare uguaglianza
% if any(abs(v_cp_y-reshape(Vcpy,1,N_tail*Nb))>1e-1) ||any(abs(v_cp_x-reshape(Vcpx,1,N_tail*Nb))>1e-1) ||any(abs(y_cp-reshape(Y_cp,1,N_tail*Nb))>1e-1) ||any(abs(x_cp-reshape(X_cp,1,N_tail*Nb))>1e-1)
% evy=any(abs(v_cp_y-reshape(Vcpy,1,N_tail*Nb))>1e-1) 
% evx=any(abs(v_cp_x-reshape(Vcpx,1,N_tail*Nb))>1e-1)
% ex=any(abs(x_cp-reshape(X_cp,1,N_tail*Nb))>1e-1)
% ey=any(abs(y_cp-reshape(Y_cp,1,N_tail*Nb))>1e-1)
% keyboard
% end
end

len = sqrt(dxp.^2 + dyp.^2);
%% 
% Normal vectors to each panel
nx = -dyp ./ len;
ny =  dxp ./ len;
vn_target = v_cp_x .* nx + v_cp_y .* ny;
% Subtracting the velocity induced by the vortices in the wake
% 

if ~isempty(vortices.pos) && t > dt
    vn_fast='si';
    switch vn_fast
        case 'no'
    for jj = 1:Nb*N_tail
        u_ind = 0;
        v_ind = 0;
        for ii = 1:size(vortices.pos,1)
            rx = x_cp(jj) - vortices.pos(ii,1);
            ry = y_cp(jj) - vortices.pos(ii,2);
            r2 =(rx^2 + ry^2);
%% 
% 
            g_sigma = 1-exp(-r2/2/sigma^2);
            u_i =  vortices.gamma(ii) * ry / (2*pi*r2)*g_sigma;
            v_i = -vortices.gamma(ii) * rx / (2*pi*r2)*g_sigma;
            u_ind = u_ind + u_i;
            v_ind = v_ind + v_i;
        end
        vn_target(jj) = vn_target(jj) - (u_ind * nx(jj) + v_ind * ny(jj));
    end
        case 'si'
% dimensioni
M = numel(x_cp);
P = size(vortices.pos,1);

% ensure column vectors
x = x_cp(:);
y = y_cp(:);
nx = nx(:);
ny = ny(:);

% create matrices (M x P): each row = control point, each col = vortex
[Xv, Xcp] = meshgrid(vortices.pos(:,1), x);   % Xv: MxP from vortices, Xcp: MxP from cp
[Yv, Ycp] = meshgrid(vortices.pos(:,2), y);

rx = Xcp - Xv;   % MxP
ry = Ycp - Yv;   % MxP

r2 = rx.^2 + ry.^2;

% avoid division by zero (points exactly at vortex position)
r2(r2 == 0) = eps;

g_sigma = 1 - exp(-r2 ./ (2 * sigma^2));

% gamma as row for broadcasting
gamma_row = vortices.gamma(:)';  % 1xP

% induced velocities from each vortex at each control point (MxP)
u_mat =  gamma_row .* (ry ./ (2*pi .* r2)) .* g_sigma;   % MxP
v_mat = -gamma_row .* (rx ./ (2*pi .* r2)) .* g_sigma;   % MxP

% sum contributions over vortices (sum across columns)
u_ind = sum(u_mat, 2);   % Mx1
v_ind = sum(v_mat, 2);   % Mx1

% update normal velocity target
vn_target = vn_target(:) - (u_ind .* nx + v_ind .* ny);
vn_target=reshape(vn_target,1,[]);
    end
end
% Computation of the influence matrix A
% The position where the trailing edge vortex forms is at one third of the distance 
% between the trailing edge of the tail and the last vortex shed.
if 0%t<=dt
    xTEV = x_tail + l(end)/Nb*cos(thetaj(end));
    yTEV = y_tail + l(end)/Nb*sin(thetaj(end));
else
    xTEV = x_tail + (vortices.pos(end,1)-x_tail)/3;
    yTEV = y_tail + (vortices.pos(end,2)-y_tail)/3;
end
%% 

% Computation of the first Nb rows of the matrix A with the relation between 
% the circulation and the normal velocity
A = zeros(Nb*N_tail+1, Nb*N_tail+1);
for ii = 1:Nb*N_tail
    for jj = 1:Nb*N_tail
        if ii == jj
            A(ii,jj) = 1/2/pi/sigma;
        else
            x_v = x_cp(jj);
           % y_v = 0.5 * (y_seg(jj) + y_seg(jj+1));
           y_v=y_cp(jj);
            rx = x_cp(ii) - x_v;
            ry = y_cp(ii) - y_v;
            r2 = sqrt(rx^2 + ry^2);
            g_sigma = 1-exp(-r2^2/2/sigma^2);
            u = -ry / (2*pi*r2^2)*g_sigma;
            v =  rx / (2*pi*r2^2)*g_sigma;
            A(ii,jj) = u * nx(ii) + v * ny(ii);
        end
    end
    rx = x_cp(ii) - xTEV;
    ry = y_cp(ii) - yTEV;
    r2 = rx^2 + ry^2;
    g_sigma = 1-exp(-r2/2/sigma^2);
    u = -ry / (2*pi*r2)*g_sigma;
    v =  rx / (2*pi*r2)*g_sigma;
    A(ii,Nb*N_tail+1) = u * nx(ii) + v * ny(ii);
end
%% 
% Computation of the last row of the matrix A imposing the Kelvin condition
A(end,:) = 1;



%% 
% Computation of the vector on the Right Hand Side of the equation
b = [vn_target'; sum(Gamma_b_old(1:end-1))];
%% 
% Solution of the equation
gamma_b = A \ b;
%% 
% Storage of the circulations
vortices.gamma_b = gamma_b;
gamma_new = gamma_b(end);
Gamma_panels = sum(gamma_b);
Gamma_vortices = sum(vortices.gamma);
gamma_total = Gamma_panels + Gamma_vortices;
% Computation of the hydrodynamic force
dGamma_dt = zeros(Nb*N_tail,1);
alpha_relax = 0.1;
if ~isempty(Gamma_b_old) && t>dt
%% 
% Derivative of the cumulative circulation from the joint to the panel " i "
% 
% 
    for jj = 1:Nb*N_tail
        dGamma_dt(jj) = (sum(gamma_b(1:jj))-sum(Gamma_b_old(1:jj)))/dt;
    end
%% 
% Relaxation of the derivative to remove spikes and improve numerical stability
    dGamma_dt = alpha_relax*dGamma_dt+(1-alpha_relax)*dGamma_dt_old';
else
    dGamma_dt = zeros(Nb*N_tail,1);
end
%% 
% Computation of the pressure difference on each panel
% 
%
nx=reshape(nx,size(len));
ny=reshape(ny,size(len));
F_vort=zeros(2,N_tail);
for I_ele=1:N_tail
Delta_p = zeros(Nb,1);
for jj = 1:Nb
    Delta_p(jj) = max(min(-rho*([v_cp_x((I_ele-1)*Nb+jj),v_cp_y((I_ele-1)*Nb+jj)]*[cos(thetaj(I_ele));sin(thetaj(I_ele))]*gamma_b((I_ele-1)*Nb+jj)/len((I_ele-1)*Nb+jj)+dGamma_dt((I_ele-1)*Nb+jj)),20000),-20000);
    if abs(Delta_p(jj)) >= 20000
        warning('t=%.3f: Unexpected pressure %.1e Pa on panel %d. Clipped at 20 kPa', t, Delta_p(jj), jj);
    end
   
end
F_vort_x = sum(Delta_p'.*len((I_ele-1)*Nb + (1:Nb)).*nx((I_ele-1)*Nb + (1:Nb)))*h;
F_vort_y = sum(Delta_p'.*len((I_ele-1)*Nb + (1:Nb)).*ny((I_ele-1)*Nb + (1:Nb)))*h;
F_vort(:,I_ele) = [F_vort_x; F_vort_y];
lj = cumsum(len((I_ele-1)*Nb + (1:Nb)));
M_vort(I_ele) = sum(Delta_p'.*len((I_ele-1)*Nb + (1:Nb)).*lj)*h;
end
%% 
% Computation of forces and moments
% 
% 
% Update of the vortex positions in the wake
% 
if ~isempty(vortices.pos) && dt > 0
    
    pos = vortices.pos;
    gamma = vortices.gamma(:);
    dx = pos(:,1) - pos(:,1).';
    dy = pos(:,2) - pos(:,2).';
    r2 = dx.^2 + dy.^2;
    g_sigma = 1-exp(-r2 / (2 * sigma^2));
    coeff = (gamma.' ./ (2 * pi * r2)) .* g_sigma;
%% 
% Avoid self induction. The coefficient on the diagonal of the coeffient matrix 
% is zero
    coeff(r2==0) = 0;
%% 
% Velocity induced by all the other vortices
    u_ind = -sum(coeff .* dy, 2);
    v_ind =  sum(coeff .* dx, 2);
%% 
% Update of the vortex positions
    new_pos = pos + dt * [u_ind, v_ind];
    vortices.pos = new_pos;
end
%% 
% Fade away the old vortices to simulate viscosity and to speed up computations
% 
% $$\Gamma_{w,k}\,(t) = \Gamma_{w,k}\,(0) \; \mathrm{e}^{-\;\frac{t}{T_{fade}}}$$
if ~isinf(T_fade)
   fade = exp(-vortices.age / T_fade);
   vortices.gamma = vortices.gamma.*fade;
end
%% 
% Store the new vortices variables
vortices.pos = [vortices.pos; xTEV, yTEV];
vortices.gamma = [vortices.gamma; gamma_new];
vortices.age = [vortices.age+dt; 0];
vortices.F_vort = [vortices.F_vort; F_vort'];
vortices.gammaTotal = [vortices.gammaTotal;gamma_total];
vortices.dGammadt = [vortices.dGammadt,dGamma_dt];
 
vort = vortices;
end