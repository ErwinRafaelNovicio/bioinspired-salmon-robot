clear all
clc 

addpath("Functions\")
addpath("Funzioni_fem\")
[A,w,p,L,area_x,a_x,MJ,K_coef]=load_paperdata;
x=(0:1:L);

E=0.03e-6;%MPa=Nmm2 , vedi grafico a p68
%rho=0.744*1e-6;%kg/mm3
rho=997*1e-9;
N=length(x);

%% material 
%      index, m,  EA ,EJ

mat=zeros(N-1,4);
mat(:,1)=1:N-1;
%mat(:,2)=1;

for I=1:N-1
    mat(I,2)=rho/(x(I+1)-x(I))*integral(@(x)area_x(x),x(I),x(I+1));%M/dL
mat(I,3)=E/(x(I+1)-x(I))*integral(@(x)area_x(x),x(I),x(I+1));%EA
    mat(I,4)=E/(x(I+1)-x(I))*integral(@(x)MJ(x),x(I),x(I+1));%EJ
end


node=[x;zeros(size(x))]';

NODE=[(1:N)',zeros(N,3),node];
ElEMENT=[(1:N-1);1:N-1;2:N;1:N-1]';


%% constraint

NODE(1,2:4)=[1,1,1];
file_name='Funzioni_fem\mesh';
input_compiler(file_name ,NODE,ElEMENT,mat)

%% lunghezza massima parte flessibile assumendo stesse proprietà corpo salmone
omegaMax=1.5;
eta=1.2;%fattore sicurezza
[EI_m_Peggiore,I_peggiore]=min(sqrt(mat(:,4)./(mat(:,2))));
Lmax=sqrt(pi^2/(eta*omegaMax)*EI_m_Peggiore);
disp(['lunghezza massima parte flessibile assumendo stesse proprietà corpo salmone:',num2str(Lmax)])