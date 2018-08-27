function references=rd_MEDLINE(fname,global_num)
% Reads a reference list from a MEDLINE (PubMed)-exported file
% tested on October 26th 2009
%
% only Journal Articles (type 1) are explicit types, all 
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
short='';
year='';
if fid~=-1,
    tline = fgetl(fid);
    nl=nl+1;
    while ischar(tline)
        [token,arg]=strtok(tline,'-');
        arg=strtrim(arg(2:end));
        tline = fgetl(fid);
        if ~ischar(tline), break; end;
        if length(tline)<4, tline='    '; end;
        while isempty(strtrim(tline(1:4))),
            arg=sprintf('%s %s',arg,strtrim(tline));
            tline = fgetl(fid);
            if ~ischar(tline); break; end;
            if length(tline)<4, tline='    '; end;
        end;
        if ~ischar(tline); break; end;
        switch strtrim(token)
            case 'PMID'
                if rec>0,
                    short=[short ':' year '_' sprintf('%i',references(rec).PMID)];
                    references(rec).short=short;
                end;
                short='';
                year='';
                in_record=1;
                rec=rec+1;
                global_num=global_num+1;
                references(rec).format=5;
                references(rec).type=1;
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
                id=str2num(arg);
                if ~isnan(id),
                    references(rec).PMID=id;
                end;
            case 'TI'
                references(rec).title=strtrim(arg);
            case 'TA'
                references(rec).journal=strtrim(arg);
            case 'DP'
                data=strtrim(arg);
                references(rec).year=data(1:4);
                year=data(1:4);
            case 'VI'
                references(rec).volume=strtrim(arg);
            case 'IP'
                references(rec).issue=strtrim(arg);
            case 'PG'
                references(rec).pages=strtrim(arg);
            case 'AU'
                if ~isempty(references(rec).authors),
                    references(rec).authors=strcat(references(rec).authors,'; ',strtrim(arg));
                else
                    references(rec).authors=strtrim(arg);
                    short=strtrim(strtok(arg));
                end;
            case 'AID'
                [doi,type]=strtok(strtrim(arg),'[');
                if strfind(lower(type),'doi'),
                    references(rec).DOI=clean_it(doi);
                end;
            case 'PT'
                if strcmpi(strtrim(arg),'Review'),
                    references(rec).type=9;
                end;
        end;
    end
    if rec>0,
        short=[short ':' year '_' sprintf('%i',references(rec).PMID)];
        references(rec).short=short;
    end;
    fclose(fid);
end;

function str=clean_it(str)
% clean up a Medline string
str=strtrim(str);
clean=0;
while ~clean,
    clean=1;
    poi=strfind(str,'&lt;');
    if ~isempty(poi),
        clean=0;
        if poi(1)<2,
            if length(str)>4,
                str=['<' str(5:end)];
            else
                str='<';
            end;
        else
            if length(str)>poi(1)+3,
                str=[str(1:poi(1)-1) '<' str(poi(1)+4:end)];
            else
                 str=[str(1:poi(1)-1) '<'];
            end;
        end;
    end;
end;
clean=0;
while ~clean,
    clean=1;
    poi=strfind(str,'&gt;');
    if ~isempty(poi),
        clean=0;
        if poi(1)<2,
            if length(str)>4,
                str=['>' str(5:end)];
            else
                str='>';
            end;
        else
            if length(str)>poi(1)+3,
                str=[str(1:poi(1)-1) '>' str(poi(1)+4:end)];
            else
                 str=[str(1:poi(1)-1) '>'];
            end;
        end;
    end;
end;
clean=0;
while ~clean,
    clean=1;
    poi=strfind(str,'lt;');
    if ~isempty(poi),
        clean=0;
        if poi(1)<2,
            if length(str)>3,
                str=['<' str(4:end)];
            else
                str='<';
            end;
        else
            if length(str)>poi(1)+2,
                str=[str(1:poi(1)-1) '<' str(poi(1)+3:end)];
            else
                 str=[str(1:poi(1)-1) '<'];
            end;
        end;
    end;
end;
clean=0;
while ~clean,
    clean=1;
    poi=strfind(str,'gt;');
    if ~isempty(poi),
        clean=0;
        if poi(1)<2,
            if length(str)>3,
                str=['>' str(4:end)];
            else
                str='>';
            end;
        else
            if length(str)>poi(1)+2,
                str=[str(1:poi(1)-1) '>' str(poi(1)+3:end)];
            else
                 str=[str(1:poi(1)-1) '>'];
            end;
        end;
    end;
end;
clean=0;
while ~clean,
    clean=1;
    poi=strfind(str,'&amp;');
    if ~isempty(poi),
        clean=0;
        if poi(1)<2,
            if length(str)>5,
                str=[str(6:end)];
            else
                str='<';
            end;
        else
            if length(str)>poi(1)+4,
                str=[str(1:poi(1)-1) str(poi(1)+5:end)];
            else
                 str=str(1:poi(1)-1);
            end;
        end;
    end;
end;
