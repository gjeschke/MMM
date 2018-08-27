function wr_rotlib_test_data(fname,testset,no_context)
% function wr_rotlib_test_data(fname,testset)
%
% reads test set for a rotamer library test
%
% fname     file name for test set
% testset   test data set structure with fields
%           pdb
% no_context optional flag that indicates context-free computation
%
% programmed for G. Jeschke, Progr. Nucl. Magn. Reson., 2013
%

if ~exist('no_context','var'),
    no_context = false;
end;

if no_context,
    cmsg = ' (no structural context considered)';
else
    cmsg = '';
end;
fid=fopen(fname,'wt');
if fid==-1,
    add_msg_board('ERROR: Output file cannot be written.');
    return;
end;

fprintf(fid,'%% Test results for library %s at %i K%s\n',testset.library,testset.T,cmsg);
fprintf(fid,'%% automatically created by MMM\n');
if testset.T~=0,
    fprintf(fid,'# TEST %s %i\n',testset.library,testset.T);
else
    fprintf(fid,'# TEST %s\n',testset.library);
end;
fprintf(fid,'%% PDB \t(ch)res1  \t(ch)res2  \t<r>/rm \tsigr \t<rsim> \tsigr_s \ttype \tsec \texp \tmode \tchainseq \tref\n');
for k=1:length(testset.data),
    md=testset.data(k);
    fprintf(fid,' %s \t%12s \t%12s \t%4.2f \t%4.2f \t%4.2f \t%4.2f \t%i \t%s \t%s \t%i',md.pdb,md.res1,md.res2,md.rexp,md.sexp,md.rsim,md.ssim,md.type,md.sec,md.exposure,md.cmode);
    if md.cmode==2,
        fprintf(fid,' \t%8s',md.chainseq);
    else
        fprintf(fid,'         ');
    end;
    fprintf(fid,' \t%s',md.ref);
    if md.cmode == 3,
        fprintf(fid,' \t\t%5.3f\t%5.3f',md.rsim_sd,md.ssim_sd);
    end;
    fprintf(fid,'\n');
end;
fprintf(fid,'# END\n\n');
fclose(fid);


% fid=fopen('rotamer_test_reduced.dat','wt');
% if fid==-1,
%     add_msg_board('ERROR: Output file cannot be written.');
%     return;
% end;
% 
% fprintf(fid,'%% Test results for library %s at %i K\n',testset.library,testset.T);
% fprintf(fid,'%% automatically created by MMM\n');
% if testset.T~=0,
%     fprintf(fid,'# TEST %s %i\n',testset.library,testset.T);
% else
%     fprintf(fid,'# TEST %s\n',testset.library);
% end;
% fprintf(fid,'%% PDB \t(ch)res1  \t(ch)res2  \t<r>/rm \tsigr \t<rsim> \tsigr_s \ttype \tsec \texp \tmode \tchainseq \tref\n');
% for k=1:length(testset.data),
%     md=testset.data(k);
%     pred = sqrt((md.uc1^2+md.uc2^2));
%     if pred <= 0.1,
%         fprintf(fid,' %s \t%12s \t%12s \t%4.2f \t%4.2f \t%4.2f \t%4.2f \t%i \t%s \t%s \t%i',md.pdb,md.res1,md.res2,md.rexp,md.sexp,md.rsim,md.ssim,md.type,md.sec,md.exposure,md.cmode);
%         if md.cmode==2,
%             fprintf(fid,' \t%8s',md.chainseq);
%         else
%             fprintf(fid,' \t        ');
%         end;
%         fprintf(fid,' \t%s\n',md.ref);
%     end;
% end;
% fprintf(fid,'# END\n\n');
% fclose(fid);
% 
fid=fopen('zhenia_tight_sites_pulse_statistics.dat','wt');
if fid==-1,
    add_msg_board('ERROR: Second output file cannot be written.');
    return;
end;
for k=1:length(testset.data),
    md=testset.data(k);
    fprintf(fid,' %4.2f \t%4.2f \t%6.4f \t%6.4f \t%6.4f\n',md.rexp,md.rsim,md.rsim_sd,md.uc1,md.uc2);
end;
fclose(fid);