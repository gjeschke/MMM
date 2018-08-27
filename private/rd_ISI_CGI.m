function references=rd_ISI_CGI(fname,global_num)
% Reads a reference list from an ISI (Web of Knowledge)-exported file
% the file must have been exported with the "Save to EndNote, RefMan, 
% ProCite" button (in *.cgi format)
% tested on October 2th 2009
%
% only Journal Articles (type 1) and Patents (8) are explicit types, all 
% other source types have type 0, but fields are still imported for them
% the short name is set to a running global reference number 
% in square brackets, which is initialized by global_num
%
% G. Jeschke, 2009

references=[];

fid=fopen(fname,'r');
nl=0;
rec=0;
in_record=0;
if fid~=-1,
    tline = fgetl(fid);
    nl=nl+1;
    while ischar(tline)
        [token,arg]=strtok(tline);
        tline = fgetl(fid);
        if ~ischar(tline), break; end;
        if length(tline)<2, tline='  '; end;
        while isempty(strtrim(tline(1:2))),
            arg=sprintf('%s; %s',strtrim(arg),strtrim(tline));
            tline = fgetl(fid);
        end;
        switch token
            case 'EF'
                break;
            case 'ER'
                pg='';
                if ~isempty(bp), 
                    pg=bp;
                elseif ~isempty(ep),
                    pg=ep;
                end;
                if ~isempty(bp) && ~isempty(ep),
                    pg=sprintf('%s-%s',bp,ep);
                end;
                references(rec).pages=pg;
                in_record=0;
            case 'PT'
                if rec>0 && isempty(references(rec).journal) && ~isempty(source),
                    references(rec).journal=source;
                end;
                source='';
                bp='';
                ep='';
                in_record=1;
                rec=rec+1;
                global_num=global_num+1;
                references(rec).format=1;
                references(rec).type=0;
                references(rec).title='';
                references(rec).book_title='';
                references(rec).authors='';
                references(rec).editors='';
                references(rec).journal='';
                references(rec).year='';
                references(rec).volume='';
                references(rec).issue='';
                references(rec).chapter='';
                references(rec).pages='';
                references(rec).DOI='';
                references(rec).publisher='';
                references(rec).city='';
                references(rec).short=sprintf('[%i]',global_num);
                references(rec).URL='';
                references(rec).PMID=0;
                type=strtrim(arg);
                switch type
                    case 'J'
                        references(rec).type=1;
                    case 'B'
                        references(rec).type=2; % book
                    case 'S'
                        references(rec).type=3; % book section
                    case 'P'
                        references(rec).type=8;
                end;
            case 'TI'
                references(rec).title=strtrim(arg);
            case 'JI'
                references(rec).journal=strtrim(arg);
            case 'SO'
                source=strtrim(arg);
            case 'SE'
                references(rec).book_title=strtrim(arg);
            case 'ED'
                references(rec).editors=strtrim(arg);
            case 'PU'
                references(rec).publisher=strtrim(arg);
            case 'AU'
                authors0=strtrim(arg);
                authors='';
                for kk=1:length(authors0),
                    if ~strcmp(authors0(kk),','),
                        authors=[authors authors0(kk)];
                    end;
                end;
                references(rec).authors=authors;
            case 'Chapter'
                references(rec).chapter=strtrim(arg);
            case 'PY'
                references(rec).year=strtrim(arg);
            case 'VL'
                references(rec).volume=strtrim(arg);
            case 'IS'
                references(rec).issue=strtrim(arg);
            case 'PN'
                all=strtrim(arg);
                if length(all)>=2,
                    references(rec).volume=all(1:2);
                end;
                if length(all)>2,
                    references(rec).issue=all(3:end);
                end;
            case 'BP'
                bp=strtrim(arg);
            case 'EP'
                ep=strtrim(arg);
            case 'DI'
                references(rec).DOI=strtrim(arg);
            case 'UR'
                references(rec).URL=strtrim(arg);
            case 'PI'
                references(rec).city=strtrim(arg);
        end;
    end
    if rec>0 && isempty(references(rec).journal) && ~isempty(source),
        references(rec).journal=source;
    end;
    fclose(fid);
end;
