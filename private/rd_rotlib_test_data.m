function testset=rd_rotlib_test_data(fname)
% function testset=rd_rotlib_test_data(fname)
%
% reads test set for a rotamer library test
%
% fname     file name for test set
% testset   test data set structure with fields
%           pdb
%
% programmed for G. Jeschke, Progr. Nucl. Magn. Reson., 2013
%

scriptfile = 'report_script.mmm';
repfid = fopen(scriptfile,'wt');
repfile = 'analysis_test_pairs.dat';

testset=[];

fid=fopen(fname);
if fid==-1,
    add_msg_board('ERROR: Test data set file does not exist.');
    return;
end;

clear testset
poi=0;
modus=0;

while 1
    tline = fgetl(fid);
    if ~ischar(tline) || modus<0, break, end
    if ~isempty(tline),
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmp(char(args(1)),'#'),
            switch upper(char(args(2)))
                case 'TEST'
                    modus=1;
                    testset.library=char(args(3));
                    if length(args)>=4,
                        testset.T=str2double(char(args(4)));
                    else
                        testset.T=0;
                    end;
                case 'REFERENCES'
                    modus=0;
                case 'END'
                    modus=-1;
                otherwise
                    modus=0;
                    fprintf(1,'Unknown keyword: %s\n',char(args(2)));
            end;
        elseif modus>0 && ~strncmp(strtrim(char(args(1))),'%',1),
            switch modus
                case 1
                    poi=poi+1;
                    testset.data(poi).pdb=char(args(1));
                    testset.data(poi).res1=char(args(2));
                    testset.data(poi).res2=char(args(3));
                    testset.data(poi).rexp=str2double(char(args(4)));
                    testset.data(poi).sexp=str2double(char(args(5)));
                    testset.data(poi).rsim=str2double(char(args(6)));
                    testset.data(poi).ssim=str2double(char(args(7)));
                    testset.data(poi).type=str2double(char(args(8)));
                    testset.data(poi).sec=char(args(9));
                    testset.data(poi).exposure=char(args(10));
                    cmode=str2double(char(args(11)));
                    testset.data(poi).cmode=cmode;
                    adr1 = sprintf('[%s]%s',char(args(1)),char(args(2)));
                    adr2 = sprintf('[%s]%s',char(args(1)),char(args(3)));
                    mk_report_script_rotlib_test(repfid,repfile,adr1,adr2)
                    if cmode==2,
                        testset.data(poi).chainseq=char(args(12));
                        testset.data(poi).ref=char(args(13));
                    else
                        testset.data(poi).ref=char(args(12));
                    end;
                    if cmode==3,
                        testset.data(poi).rsim_sd = 0;
                        testset.data(poi).ssim_sd =0;
                        testset.data(poi).uc1 =0;
                        testset.data(poi).uc2 =0;
                    end;
            end;
        end;
    end;
end;

fclose(fid);
fclose(repfid);
