clc
close all
clear all


%dati
a=[21,21,16,11.5];%[mm] distanza fori centro vertebre
b=[70,25,25];%[mm]
c=[25,25,25];%[mm]

r_scotch_yoke=8*1e-3;
%
anima=true;

%______________________________________________
rot=@(a) [cos(a),-sin(a);sin(a),cos(a)];
l0=zeros(3,1);
%
%plot movimento di sola mezza coda
figure;
col={'b','r','g','y'};
axis([-0.1,sum(b+c),-5*a(1),5*a(1)])
for I=1:3
[l0(I),lam(I)]=initial_len(a(I),a(I+1),b(I),c(I));
h(I) = animatedline('LineWidth', 2,'Color',col{I});
end
h(4)= animatedline('LineWidth', 2,'Color',col{4},'LineStyle','none','marker','o');
omega=2*pi*2;
t=0:0.1/omega:3;
alpha_k=zeros(4,length(t));





%% imposta ampiezze

T_fish=[   0.0249    0.0287    0.0253]';%ampiezze coppie motrici prese da dvm
height=25; % spostamento verticale ultima vertebra attuata
s=@(dl) spos2(dl,a,b,c,l0);
dl_star=fzero(@(dl) height-s(dl),0.2)
[y,Amp_vec]=spostamento(dl_star,a,b,c,l0)

%% animazione 
yline(height);

if anima
for K=1:length(t)

dl=dl_star*sin(omega*t(K));%dl function
store=zeros(2,7,3);
alpha=zeros(4,1);

for I=1:3
[alpha_suc,points]=base(a(I),a(I+1),b(I),c(I),l0(I)-dl*I,'animazione');


alpha(I+1) = alpha_suc; % Store the computed angle for each iteration

trans=eye(2);
for j=1:I
trans=trans*rot(alpha(j));
end
points=trans*points;

if I>1
store(:,:,I)=store(:,3,I-1)+points;
else
    store(:,:,I)=points;
end

end
alpha_k(:,K)=alpha;

%%

hold on
xlabel('X-axis')
ylabel('Y-axis')
title('Plot of Points')
grid on

for I=1:3+1
clearpoints(h(I));
end
for I=1:3
     %plot([store(1,:,I),store(1,1,I)],[store(2,:,I),store(2,1,I)])
 addpoints(h(4),store(1,2,I),store(2,2,I));
 addpoints(h(I), [store(1,:,I),store(1,3,I)], [store(2,:,I),store(2,3,I)]);
end

drawnow
pause(.01)
end
alpha_n=alpha_k./(alpha_k(2,:)+eps);
end

figure
subplot(2,1,1)
title('rel angle')
plot(t,(alpha_k))
subplot(2,1,2)
title('abs angle')
plot(t,cumsum(alpha_k))
%% lagrange
save_dlam=zeros(3,length(t));
save_dalpha=zeros(3,length(t));
for k=1:length(t)
    dl=dl_star*sin(omega*t(k));
for I=1:3
Delta_li=l0(I)+dl*I;
[dalpha(I),~,dlam(I)]= base(a(I),a(I+1),b(I),c(I),Delta_li);
d=sqrt(a(I)^2+b(I)^2);
e=sqrt(a(I+1)^2+c(I)^2);
cost(I)=-2*dl/sqrt(1-((d^2+e^2-Delta_li^2)/(2*d*e))^2);
end
dalpha=dalpha+0*pi/9*[0,1,2];%sfasamento
save_dalpha(:,k)=dalpha;
T1=1/r_scotch_yoke*sin(omega*t(k)); %tiro fune per unità di coppia motore

%Tqi=  T_passante+T_termina_nella_vertebra=Cm/r*sin(wt)*d/dalpahai(punto
%applicazione) 
Tp1=2*sum([1/2,1/3])*T1*[0;cos(dlam(1))*tan(dlam(2))-sin(dlam(1))]'*rot(dalpha(1))*rot(pi/2)*[c(1);a(2)]*1e-3;
Ttnv1=T1*[-cos(dlam(1));cos(dlam(1))*tan(dlam(2))-sin(dlam(1))]'*rot(dalpha(1))*rot(pi/2)*[c(1);a(2)]*1e-3;


%
Tp2=1/3*T1*cos(dlam(1))/cos(dlam(2))*[0;cos(dlam(2))*tan(dlam(3))-sin(dlam(2))]'*(rot(dalpha(2))*rot(pi/2)*[c(2);a(3)])*1e-3;
Ttnv2=T1*[-cos(dlam(1));cos(dlam(1))*tan(dlam(3))-sin(dlam(2))]'*(rot(dalpha(2))*rot(pi/2)*[c(2);a(3)])*1e-3;

%
Tp3=T1*[0;cos(dlam(1))*tan(dlam(3))]'*(rot(dalpha(3))*rot(pi/2)*[c(3);a(4)])*1e-3;
Ttnv3=0;




 %Tqi=T1*(dpos/dalphai)      
Tq1(k)=Tp1+Ttnv1;
Tq2(k)=Tp2+Ttnv2;
Tq3(k)=Tp3+Ttnv3;
   

end

tiro_lagrange=T_fish'./rms([Tq1;Tq2;Tq3]');
figure
subplot(2,1,1)
hold on
plot(t,tiro_lagrange(1)*Tq1)
plot(t,tiro_lagrange(2)*Tq2)
plot(t,tiro_lagrange(3)*Tq3)
legend('Tq1','Tq2','Tq3')
subplot(2,1,2)
hold on
plot(t,save_dalpha)
legend('alfa1','alfa2','alfa3')



Cm_lag=mean(T_fish'./rms([Tq1;Tq2;Tq3]'))
Pm_lag=Cm_lag*omega
%% functions
function [alpha_suc,points,lam]= base(a1,a2,b,c,dl,caso)
if nargin < 6 || isempty(caso)
    caso = '';
end
d=sqrt(a1^2+b^2);
e=sqrt(a2^2+c^2);
gamma=atan(a1/b);
delta=atan(a2/c);
rho=acos((d^2+e^2-dl^2)/(2*d*e));

alpha_suc=pi-(delta+gamma+rho);
rot=@(a) [cos(a),-sin(a);sin(a),cos(a)];
points=zeros(2,5);
points(:,1)=[0;0];%A,B,c,D,E
points(:,2)=[b;0];
points(:,3)=[b;0]+rot(alpha_suc)*[c;0];
points(:,4)=[b;0]+rot(alpha_suc+delta)*[e;0];
points(:,5)=[0,a1];
if strcmp(caso, 'animazione')
points(:,7)=[b;0]+rot(alpha_suc-delta)*[e;0];%D'
points(:,6)=[0,-a1];%E'
end
lam=acos((dl^2+e^2-d^2)/(2*dl*e))-delta;
end

%%
function [l0,lam]= initial_len(a1,a2,b,c)
d=sqrt(a1^2+b^2);
e=sqrt(a2^2+c^2);
gamma=atan(a1/b);
delta=atan(a2/c);
l0=sqrt(d^2+e^2-2*e*d*cos(pi-(gamma+delta)));
lam=acos((l0^2+e^2-d^2)/(2*l0*e))-delta;
end
%%
function [y,alpha]=spostamento(dl,a,b,c,l0)
rot=@(a) [cos(a),-sin(a);sin(a),cos(a)];
store=zeros(2,5,3);
alpha=zeros(4,1);

for I=1:3
[alpha_suc,points]=base(a(I),a(I+1),b(I),c(I),l0(I)-I*dl);


alpha(I+1) = alpha_suc; % Store the computed angle for each iteration
trans=eye(2);

for j=1:I
trans=trans*rot(alpha(j));

end
points=trans*points;

if I>1
store(:,:,I)=store(:,3,I-1)+points;

else
    store(:,:,I)=points;
end

end
y=store(2,4,end);%prendo C
%y=store(2,3,end);%prendo B
alpha=alpha';
end
function y=spos2(dl,a,b,c,l0)


[y,~]=spostamento(dl,a,b,c,l0);
end
