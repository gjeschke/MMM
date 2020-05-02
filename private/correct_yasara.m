function message = correct_yasara(fname,substitute,newcid,oname)

if ~exist('substitute','var')
    substitute = '';
end
if ~exist('newcid','var')
    newcid = '';
end


message.error = 0;
message.text = 'OK';

ifid=fopen(strcat(fname,'.pdb'));
if ifid==-1
    message.error=1;
    message.text='Input file could not be opened';
    return;
end

if exist('oname','var') && ~isempty(oname)
    ofid=fopen(oname,'wt');
else
    ofid=fopen(strcat(fname,'_corr.pdb'),'wt');
end
if ofid==-1
    message.error=2;
    message.text='Output file could not be written';
    return;
end

chain = '';
tline = '';
while 1
    tline0 = tline;
    tline = fgetl(ifid);
    if ~ischar(tline), break, end
    if length(tline)<16 % catches too short end line of MolProbity files
        fprintf(ofid,'%s\n',tline);
    else
        record = tline(1:6);
        % fprintf(1,'%s\n',record);
        if strcmp(record,'HETATM')
            record = 'ATOM  ';
        end
        if strcmp(record,'ATOM  ')
            for ks = 1:length(substitute)
                if strcmpi(tline(22),substitute(ks))
                    tline(22) = upper(newcid(ks));
                end
            end
            if ~isempty(chain)
                if ~strcmpi(chain,tline(22))
                    tline0(1:6) = 'TER   ';
                    fprintf(ofid,'%s\n',tline0(1:27));
                end
            end
            chain = tline(22);
            atname = tline(13:16);
            skip = false;
            switch atname
                case '1H5*' 
                    atname = ' H5*';
                case '2H5*' 
                    atname = 'H5**';
                case '*HO2' 
                    atname = 'HO2*';
                case ' O1P' 
                    atname = ' OP1';
                case ' O2P' 
                    atname = ' OP2';
                case '*HO3' 
                    skip = true;
                case '*HO5' 
                    skip = true;
            end
            if ~skip
                for k = 1:length(atname)
                    if atname(k) == '*'
                        atname(k) = '''';
                    end
                end
                tline(1:6) = record;
                tline(13:16) = atname;
                fprintf(ofid,'%s\n',tline);
            end
        end
    end
end

fclose(ifid);
fclose(ofid);