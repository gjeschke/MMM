function mk_Tinker_key(network,xyzname,options,correspondence)
% function mk_Tinker_key(network,xyzname,correspondence,options)
%
% updates a Tinker key file with position constraints from the current 
% C_alpha network
%
% network   array of Calpha coordinates of network points
% xyzname   xyz file name
% options   key file options, all fields optional
%           .solvation  Tinker solvation keyword
% correspondence    correspondence table for xyz file atom numbers to MMM
%                   internal variables, prepared with assign_and_key.m
%                   .network  Tinker indices of C alpha network points
%                   .chains   array of Tinker index vectors for all chains
%                             in structure snum, .chains(kc).pointers
%                   .model    model number in structure snum, is either km
%                             or 1 (default)
%
% G. Jeschke, 2012

[mypath,basname] = fileparts(xyzname);

keyfile=fullfile(mypath,sprintf('%s.key',basname));
restrained=correspondence.network;
% write key file
fid=fopen(keyfile,'wt');

if isfield(options,'solvation') && ~isempty(options.solvation),
    fprintf(fid,'SOLVATE %s\n',options.solvation);
end;
if ~isfield(options,'restraint_force'),
    options.restraint_force=100; % this is the Tinker default
end;
[mr,~]=size(restrained);
for k=1:mr,
    fprintf(fid,'RESTRAIN-POSITION %i %12.6f%12.6f%12.6f%10.2f\n',restrained(k),network(k,1),network(k,2),network(k,3),options.restraint_force);
end;

fclose(fid);

add_msg_board(sprintf('Key file written: %s',keyfile));
