clc
close all
clear all

addpath('.\Functions\')

%% dal paper
[A,w,p,L,area_x,a_x,MJ_x,K_coef,b_x]=load_paperdata;
rho=0.744*1e-6; %kg/mm3
step0=11;
L0=linspace(0,L,step0);%[mm]
w_mean=[mean(w(:,1)),mean(w(:,2)),mean(w(:,3))];%[rad/Hz]
%% animazione
do_animation=0 ;%1 per animare
t=0:0.01:10;
steps=300;
Lv=linspace(0,L,steps);%[mm]
Lq=round(linspace(1,steps-1,10));
animation_speed=1;%velocità animazione
magnification=10;% aumento ampiezze
%% design

xgrid_ideal=[0,0.1,0.20,0.5,0.6,0.7,0.8,1]*L;
xgrid_index= round(xgrid_ideal/mean(diff(Lv)))+1;

pinna_flessibile='flessibile';% mettere 'rigida' o 'flessibile'
% determino i punti effettivi dove vado a determinare i moduli della
% struttura, dipende dal valore che do' a "steps"
switch pinna_flessibile
    case 'flessibile'
discret_coda=10;
last_rigid_el=numel(xgrid_index)-1;
xgrid_index=[xgrid_index(1:end-1),xgrid_index(end-1)+1:round(diff(xgrid_index(end-1:end))/discret_coda):xgrid_index(end)-1,xgrid_index(end)];
xgrid=Lv(xgrid_index);% posizione effettiva delle sezioni, dipende da quanti elementi ha Lv, maggiori saranno e più vicino sono ai valori di xgrid che scelgo 
   
    case 'rigida'
last_rigid_el=numel(xgrid_index);
xgrid=Lv(xgrid_index); % posizione effettiva delle sezioni, dipende da quanti elementi ha Lv, maggiori saranno e più vicino sono ai valori di xgrid che scelgo 
end

x_rigid=xgrid(1:last_rigid_el);%hinges coordinates of rigid element

%%
y=zeros(steps,numel(t));
Y_3f=zeros([size(y),3]);
for J=1:3
    %phase=@(x) interp1(L0,p(:,J),x,"linear","extrap");
    phase=@(x) polyval(polyfit(L0,p(:,1),3),x);
    Ax=@(x) polyval(polyfit(L0,A(:,J),2),x);

    figure
    subplot(2,1,1)
    plot(Lv,Ax(Lv))
    title('ampiezze')
subplot(2,1,2)
    plot(Lv,phase(Lv))
title('phase as function of x')
%%

y=zeros(steps,numel(t));

for I=1:length(t)
t4=t(I);
    y(:,I)=Ax(Lv).*sin(w_mean(J)*t4-phase(Lv));
end
Y_3f(:,:,J)= y;
if do_animation==1
    Nome=['GIF_speed',num2str(J),'.gif'];
 animazione_quadrati(Lv, magnification*y,xgrid,xgrid_index,L/steps, numel(t),animation_speed,Nome)% notare che sto  ingrandendo l'ampiezza del moto 
 sgtitle(['magnification',num2str(magnification)])
end
end

%% mass - moment inertia - c.o.m.

[X_g,M]=com(area_x,a_x,xgrid(1:last_rigid_el),rho);
MJ=moment_inertia(MJ_x,xgrid(1:last_rigid_el));%mm4
[~,M_all]=com(area_x,a_x,xgrid,rho);
Mtot=sum(M_all);
%Mtot=rho*integral(area_x,0,L);
Mrel=M/Mtot




%% molle concentrate

theta=zeros([last_rigid_el,numel(t),3]);
for J=1:3
y_model_rigid=Y_3f(xgrid_index(1:last_rigid_el),:,J);
N_rigid=numel(xgrid_index(1:last_rigid_el));

figure
for I=1:N_rigid
subplot(N_rigid,1,I)
if I==1
    a=1;%derivata in avanti
else 
    a=-1;%derivata all'indietro
end
theta(I,:,J)=atan2((y_model_rigid(I,:)+a*y_model_rigid(I+a,:)),(xgrid(I)+a*xgrid(I+a)+eps));
    plot(t,theta(I,:,J))
end
sgtitle('theta parti rigide')
end

%% stiffness

Lrobot=500;%[mm]
ratioL=Lrobot/L;
L_fea=  [x_rigid(1:end-1),linspace(x_rigid(end),L,5+1)];
kFea2=FEA3(L_fea(1:end),1e-6,7);
k_DVM=kFea2(2:end);
figure
bar(L_fea(2:end)/max(L_fea),k_DVM)
title('stifness')

%% DVM 

amm=@(x)1e3*a_x(x);%[mm]
bmm=@(x)1e3*b_x(x);%[mm]

for I=1:length(L_fea)-1
%Area_DVM(I)=integral(@(x)  pi/2*(3*(amm(x)+bmm(x)) - sqrt((3*amm(x) + bmm(x)).*(amm(x) + 3*bmm(x)))),L_fea(I),L_fea(I+1)) *ratioL^2*1e-6;%[m^2]
Area_DVM(I)=integral(@(x)  amm(x),L_fea(I),L_fea(I+1)) *ratioL^2*1e-6;%[m^2]
end
h_DVM=1/L*integral(@(x) bmm(x),0,L)*1e-3;%[m]
[~,M_DVM]=com(area_x,a_x,L_fea,997*1e-9);
M_DVM=M_DVM/sum(M_DVM)

L_DVM=L_fea*ratioL*1e-3; %[m]
dL_DVM=diff(L_DVM);
disp(['N_link: ',num2str(numel(dL_DVM))])
disp(['Lengths [m]: ',num2str(dL_DVM)])
disp(['area [m^2]: ',num2str(Area_DVM)])
disp(['h fin and body thickness [m]: ',num2str(h_DVM)])
disp(['k [Nm/rad]: ',num2str(k_DVM)])
