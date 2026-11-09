function a=input_compiler(file_name,nodes,element,mat)
string=[file_name,'.inp'];
fid=fopen(string,"w");


%% file name
fprintf(fid,"! %s\n",file_name);
%% node writing
fprintf(fid,"! ----------------------------------------------\n");
fprintf(fid,"! list of nodes:\n");
fprintf(fid,"*NODES\n");
fprintf(fid,"! n of nodes-constraint (x,y,phi)- Xnode Ynode\n");
fprintf(fid,'%d\t %d %d %d\t %f %f\n',nodes');
fprintf(fid,"*ENDNODES\n");

%% element writing
fprintf(fid,"! ----------------------------------------------\n");
fprintf(fid,"! list of elements\n");
fprintf(fid,"*BEAMS\n");
fprintf(fid,"! elem,node1,nodo2,mat\n");
fprintf(fid,"%d %d %d %d\n",element');
fprintf(fid,"*ENDBEAMS\n");


%% materials
fprintf(fid,"! ----------------------------------------------\n");
fprintf(fid,"! list of material\n");
fprintf(fid,"*PROPERTIES\n");
fprintf(fid,"! N,n,EA,EJ\n");
fprintf(fid,"%d %f %f %f\n",mat');
fprintf(fid,"*ENDPROPERTIES\n");


fclose(fid);
open(string)
a='done';
end
