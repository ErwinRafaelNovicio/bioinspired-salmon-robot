function [k,dtheta_reldq]=FEA3(x_rigid,F0,index_force)
input_file='Funzioni_fem\mesh';
addpath("Funzioni_fem\")
%% load structure
[file_i,xy,nnod,sizee,idb,ndof,incid,l,gamma,m,EA,EJ,posiz,nbeam,pr]=loadstructure_mod(input_file);

%% draw structure per vedere se sono giusti vincoli e mesh
dis_stru(posiz,l,gamma,xy,pr,idb,ndof);


%% assembly
[M,K] = assem(incid,l,m,EA,EJ,gamma,idb);%costruisco matrici massa e rigidezza
Kfem=K(1:ndof,1:ndof);%elimino gradi di libertà vincolati
dx=mean(diff(xy(:,1)));

% definisco dove applicare la forzante
iF=zeros(size(index_force));
for I=1:length(index_force)
%iF=[find(abs(xy(:,1)-x_rigid(5))<dx,1),find(abs(xy(:,1)-x_rigid(6))<dx,1),find(abs(xy(:,1)-x_rigid(7))<dx,1)];%dove applicherò la forza sul pesce
iF(I)=find(abs(xy(:,1)-x_rigid(index_force(I)))<dx,1);

end
F=zeros(ndof,1);
F(idb(iF,3))=F0;

%calcolo deformazioni
%def=Kfem\F;
def=solveFEM(Kfem,F);
y=def(2:3:end);
theta=def(3:3:end);

figure
plot([0;xy(2:end,1)],[0;y],LineWidth=2)
title('Deformed beam')

theta_x=@(x) interp1([0;xy(2:end,1)],[0;theta],x,"linear","extrap");
theta=interp1([0;xy(2:end,1)],[0;theta],x_rigid,"linear","extrap");

Y_x=@(x) interp1([0;xy(2:end,1)],[0;y],x,"linear","extrap");
Y=interp1(xy(2:end,1),y,x_rigid(1:end-1),"linear","extrap");

% load applied to robot
F_robot=zeros(size(x_rigid));
F_robot(index_force)=F0;

%dtheta_reldq, per ottenere K dei corpi rigidi;
dtheta_reldq=zeros(length(x_rigid));
dtheta_reldq(1,1)=1;
for I=2:length(x_rigid)
    dtheta_reldq(I,(1:2)+I-2)=[1,-1];
end

%d'diag(k)*d*theta=F-> theta'k'=f',
dL=diff([x_rigid(1),x_rigid])';

%k0=1e-5*ones(size(x_rigid));%Guess iniziale
k_flex=15e-4;
 k0=[1,1.9, 1.52, 0.47, 0.34, 0.20, 0.15, 10*k_flex, 5*k_flex, 2.5*k_flex,1.1* k_flex, k_flex]; 

 %minimizzo errore tra Ypesce-Ymodello e THETApesce-THETAmodello
k=fmincon(@(k) fun_error(k,theta_x,Y_x,dL,F_robot',dtheta_reldq),k0);
k=max(k,1e-6);
 K = dtheta_reldq' * diag([k']) * dtheta_reldq;   % ricostruisco K 
 
 
 residuo = norm(K*theta' - F_robot);% verifico che alla fine errore su theta sia piccolo
 disp('residuo ||K*theta - F|| ='); disp(residuo)

 %plot risultati
hold on

thetaR=K\F_robot';
Yres=[cumsum(sin(thetaR).*dL)];
Xres=[cumsum(dL.*cos(thetaR))];


plot(Xres,Yres,LineWidth=2)
plot(Xres,Yres,'o',LineWidth=2)
plot(Xres,Y_x(Xres),'square',LineWidth=2)
legend('Salmon fish model','Robot link','Robot joints')

figure
title('angle variation')
hold on 
plot(xy(2:end,1),def(3:3:end),LineWidth=2)
plot(Xres,thetaR,LineWidth=2)
plot(Xres,thetaR,'o',LineWidth=2)
plot(Xres,theta_x(Xres),'square',LineWidth=2)

legend('Salmon fish model','Robot link','Robot joints')

% p(k)
 end

function err=fun_error(k,theta,Y,dL,F,d)
k=max(k,1e-6);
 K = d' * diag([k']) * d ;  % ricostruisco K per controllo
%K
a=8;
Wy=diag([40*ones(size(1:2))';50*ones(size((2+1):a))';1*ones([numel(k)-a,1])]);
Wy(1,1)=1;

Wtheta=diag([20*ones(size(1:2))';40*ones(size((2+1):a))';1*ones([numel(k)-a,1])]);%less interested in getting the right angle for the flexible bit
Wtheta(1,1)=1;
thetaR=K\F;
%thetaR=solveFEM(K,F);

X=cumsum(cos(thetaR).*dL);
Yr=cumsum(sin(thetaR).*dL);

Y=Y(X);
theta=theta(X);

err=(theta-thetaR)'*Wtheta*(theta-thetaR)+(Y-Yr)'*Wy*(Y-Yr);

end