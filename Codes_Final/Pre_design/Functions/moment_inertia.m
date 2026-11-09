function J=moment_inertia(MJ_x,xrigid)

N_rigid=numel(xrigid)-1;
dx=diff(xrigid);
J=zeros(N_rigid,1);
for I=1:N_rigid

J(I)=integral(MJ_x,xrigid(I),xrigid(I+1))/dx(I);%media integrale

end