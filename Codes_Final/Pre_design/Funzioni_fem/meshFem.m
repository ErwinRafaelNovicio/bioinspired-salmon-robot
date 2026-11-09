function   [node,ID]=meshFem(Xdis,dL,ID_start,keep)
%INPUT 
% -Xdis position of geometrical discontinuity and point of application of force
% -dL step size
% -ID_Start first index of Xdis(1)  
% *keep if set to zero doesn't include Xdis points in nodes
%OUTPUT
%-node position of nodes
%-ID position of Xdis
if keep==0
node=[];
else
node(1,:)=Xdis(1,:);
end
  ID(1)=ID_start;
for I=1:length(Xdis)-1
    temp=[];
 
    L=norm(Xdis(I,:)-Xdis(I+1,:),2);
    dL=min(dL,L);
    phi=angle(Xdis(I+1,1)-Xdis(I,1)+1i*(Xdis(I+1,2)-Xdis(I,2)));
    if (Xdis(end,:)-Xdis(1,:))*[cos(phi);sin(phi)]<0
        phi=phi+pi;
    end
    versor=dL*[cos(phi),sin(phi)];
    N=1;    
    temp(N,:)=Xdis(I,:)+versor;
    dist=norm(temp(end,:)-Xdis(I+1,:),2);
    while dist>dL
      N=N+1;
      temp(N,:)=temp(N-1,:)+versor;
      dist=norm(temp(end,:)-Xdis(I+1,:),2);
    end
    %% avoid point duplication
    if sum( temp(end,:)==Xdis(I+1,:))>1
       if min(size(temp))==1
           temp=[];
       else
       temp=temp(1:end-1,:);
       end
       N=N-1;
    end
    if keep==0
    node=[node;temp];
    ID(I+1)=ID(I)+N;
    else
    node=[node;temp;Xdis(I+1,:)];
    ID(I+1)=N+1+ID(I);
end
end
end