function [Xg,M]=com(area,a,xrigid,rho)
%% assumo sezioni simmetriche lungo z,y ma non lungo x

N_rigid=numel(xrigid)-1;
Xg=zeros(N_rigid-1,1);
dx=diff(xrigid);
m1=@(x) x.*area(x);

for I=1:N_rigid
   V(I)=integral(area,xrigid(I),xrigid(I+1));
    
     Xg(I)=1/V(I)*integral(m1,xrigid(I),xrigid(I+1));
  
end
M=V.*dx*rho;
figure

plot(xrigid,0*ones(size(xrigid)))
hold on;
plot(xrigid,0*ones(size(xrigid)),'o')
plot(xrigid,[a(xrigid);-a(xrigid)],'r--')
xline(Xg,'--')
title('C.o.Ms rigid element ')
end