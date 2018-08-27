function wr_rotamer_statistics(fid,adr,calib_T,rotamer_stats,threshold,T,forgive,label)
% function wr_rotamer_statistics(fid,rotamer_stats,threshold,T,forgive,label)
%
% writes rotamer statistics to an ASCII file with file identifier fid
% a record for one residue with MMM address adr is written
% 
% further input parameters:
% calib_T           the calibration temperature of the rotamer library
% rotamer_stats     the rotamer statistics as computed by get_rotamers.m
% threshold         the threshold for total population of all neglected
%                   rotamers
% T                 labeling temperature
% forgive           "forgive factor"
% label             label name, is translated to a label type number (LTN)
%                   R1A (MTSL): 1
%                   IA1 (IA-PROXYL): 2
%
% the record format is:
%
% the first line is formally a Matlab comment line
% line 1:       % address calibration temperature of library label
% the following lines are ASCII data lines with five numbers per line, the
% numbers are separated by space, not commas
% line 2:       number of rotamers (NR), number of significant rotamers, T,
%               forgive, label type (LTN)
% line 3:       x, y, z coordinate of C-alpha of the residue, external net
%               potential, partition function
% line 4..3+NR: x, y, z coordinate of N-O midpoint, relative populations 
%               due to external potential, populations due to total potential

LTN=0;
switch label
    case 'R1A'
        LTN=1;
    case 'IA1'
        LTN=2;
end;
ext_pop=rotamer_stats.all_potentials.ext_pop;
nrot=length(ext_pop);
NOall=rotamer_stats.NOall;
CA=rotamer_stats.loc_frame_Ca;
cumul=cumsum(sort(rotamer_stats.all_potentials.pop_rot,'descend'));
signum=length(cumul(cumul<1-threshold))+1;
fprintf(fid,'%% %s labeled at calibration temperature %5.1f (label %s)\n',adr,calib_T,label);
fprintf(fid,'%8i%8i%8.1f%8.2f%8i\n',nrot,signum,T,forgive,LTN);
fprintf(fid,'%8.3f%8.3f%8.3f%8.3f%8.5f\n',CA(1),CA(2),CA(3),rotamer_stats.all_potentials.ext_net,rotamer_stats.all_potentials.partition_function);


% % sort NOall according to the total rotamer potential (NOall(:,4) column)--
% keyboard
NOall=sortrows(NOall,-4);

for k=1:nrot,
    fprintf(fid,'%8.3f%8.3f%8.3f%8.3f%8.5f%8.0f\n',NOall(k,1),NOall(k,2),NOall(k,3),ext_pop(k),NOall(k,4),NOall(k,5));
end;
%--------------------------------------------------------------------------

% % original Gunnar's version (07122010)
% for k=1:nrot,
%     fprintf(fid,'%8.3f%8.3f%8.3f%8.3f%8.5f\n',NOall(k,1),NOall(k,2),NOall(k,3),ext_pop(k),NOall(k,4));
% end;
