function [centrali, cerniera_primo, cerniera_secondo] = vincoli_elementi(input_file)
addpath('.\Input\')
% input_file='carlo';

eval(['fid_i=fopen(''',input_file,'.inp'',''r'');']);

findcard(fid_i,'*NODES')
finewhile=1;
iconta=0;
while finewhile ==1
    line=scom(fid_i);
    finewhile=isempty(findstr(line,'*ENDNODES'));
    if finewhile == 1
        tmp=sscanf(line,'%i %i %i %i %f %f')';
        iconta=iconta+1;
        if iconta ~=tmp(1)
        disp('Errore: nodi non numerati in ordine progressivo')
        break  
        end
        ivinc(iconta,:)=tmp(2:4);
        xy(iconta,:)=tmp(5:6);
    end
end

cerniere=find(xy(:,2)==0);     % (numero riga = numero nodo)   

findcard(fid_i,'*BEAMS')
finewhile=1;
iconta=0;
while finewhile ==1
    line=scom(fid_i);
    finewhile=isempty(findstr(line,'*ENDBEAMS'));
    if finewhile == 1
        tmp=sscanf(line,'%i %i %i %i')';
        iconta=iconta+1;
        if iconta ~=tmp(1)
        disp('Errore: nodi non numerati in ordine progressivo')
        break  
        end
        nodi(iconta,:)=tmp(2:3);
    end
end
nbeam=iconta;

centrali=[];
cerniera_primo=[];
cerniera_secondo=[];

cerniera_primo=[cerniera_primo,(find(nodi(:,1)==cerniere(1)))'];    %(le cerniere sono 2 (prima deve essere in 0,0), faccio così altrimenti prende anche nuovi nodi a y=0 ma elementi nuovi messi dopo nell' input file)
cerniera_primo=[cerniera_primo,(find(nodi(:,1)==cerniere(2)))'];

cerniera_secondo=[cerniera_secondo,(find(nodi(:,2)==cerniere(1)))'];
cerniera_secondo=[cerniera_secondo,(find(nodi(:,2)==cerniere(2)))'];

centrali=1:nbeam;
cerniere_ingiro=[cerniera_primo,cerniera_secondo];
centrali(cerniere_ingiro)=[];


end 








