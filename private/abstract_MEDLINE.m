function outname=rd_MEDLINE(fname,global_num)
% Reads references from a MEDLINE (PubMed)-exported file
% and reformats them as a HTML file for abstract display
% tested on November 11th 2009
%
%
% G. Jeschke, 2009

global general

outname=strcat(general.tmp_files,'PubMed_abstract.html');

fout=fopen(outname,'w');
fprintf(fout,'<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"\n');
fprintf(fout,'"http://www.w3.org/TR/html4/loose.dtd">\n<html>\n');

fid=fopen(fname,'r');
nl=0;
rec=0;
auon=0;
in_record=0;
authors='';
abstract='';
title='No title available';
journal='unknown journal';
volume='';
year='';
pages='';
issue='';
keywords='';
doi='';
is_review=0;
date='';
address='';
language='unknown';
references='?';
if fid~=-1,
    tline = fgetl(fid);
    nl=nl+1;
    while ischar(tline)
        if strfind(tline,'-'),
            abon=0;
        elseif abon,
            abstract=[abstract ' ' strtrim(tline)];
        end;
        [token,arg]=strtok(tline,'-');
        arg=strtrim(arg(2:end));
        tline = fgetl(fid);
        if ~ischar(tline), break; end;
        if length(tline)<4, tline='    '; end;
        while isempty(strtrim(tline(1:4))),
            arg=sprintf('%s %s',arg,strtrim(tline));
            tline = fgetl(fid);
        end;
        switch strtrim(token)
            case 'PMID'
                auon=0;
                id=str2num(arg);
                if ~isnan(id),
                    fprintf(fout,'<head>\n<title>PubMed record %i</title>\n</head>\n',id);
                else
                    fprintf(fout,'<head>\n<title>PubMed record with unknown ID</title>\n</head>\n');
                end;
                in_record=1;
                fprintf(fout,'<body>\n');
            case 'TI'
                auon=0;
                title=strtrim(arg);
            case 'TA'
                auon=0;
                journal=strtrim(arg);
            case 'VI'
                auon=0;
                volume=strtrim(arg);
            case 'IP'
                auon=0;
                issue=strtrim(arg);
            case 'PG'
                auon=0;
                pages=strtrim(arg);
            case 'DP'
                auon=0;
                data=strtrim(arg);
                year=data(1:4);
            case 'DP'
                date=strtrim(arg);
            case 'AD'
                address=strtrim(arg);
            case 'LA',
                switch lower(strtrim(arg))
                    case 'eng'
                        language='English';
                    case 'fre'
                        language='French';
                    case 'jpn'
                        language='Japanese';
                    case 'rus'
                        language='Russian';
                    case 'ger'
                        language='German';
                    case 'spa'
                        language='Spanish';
                    case 'chi'
                        language='Chinese';
                    case 'pol'
                        language='Polish';
                    case 'por'
                        language='Portugese';
                    otherwise
                        language=strtrim(arg);
                end;
            case 'AU'
                auon=1;
                if ~isempty(authors),
                    authors=[authors,'; ',strtrim(arg)];
                else
                    authors=strtrim(arg);
                end;
            case 'AID'
                auon=0;
                [doi,type]=strtok(strtrim(arg),'[');
                if strfind(lower(type),'doi'),
                    doi=strtrim(doi);
                end;
            case 'AB'
                abstract=strtrim(arg);
                abon=1;
            case 'RF'
                references=strtrim(arg);
            case 'MH'
                if isempty(keywords),
                    keywords=strtrim(arg);
                else
                    keywords=strcat(keywords,', ',strtrim(arg));
                end;
            case 'PT'
                if strcmpi(strtrim(arg),'Review'),
                    is_review=1;
                end;

        end;
    end
    fclose(fid);
else
    outname='';
end;

if isempty(authors), authors='anonymous'; end;

fprintf(fout,'<h3>%s</h3>',title);

fprintf(fout,'<p>%s <b>%s</b>, <i>%s</i>: %s(%s), %s.</p>\n',authors,year,journal,volume,issue,pages);
if ~isempty(date),
    fprintf(fout,'<p><i>Publication date: </i>%s</p>\n',date);
end;
if ~isempty(address),
    fprintf(fout,'<p>%s</p>\n',address);
end;
if ~isempty(doi),
    fprintf(fout,'<p><i>DOI</i>: %s</p>\n',doi);
end;
if isempty(abstract),
    fprintf(fout,'<p><b>No abstract available.</b></p>\n');
else
    fprintf(fout,'<h4>Abstract</h4>\n');
    fprintf(fout,'<p>%s</p>\n',abstract);
end;
if ~isempty(keywords),
    fprintf(fout,'<h4>Keywords</h4>\n');
    fprintf(fout,'<p>%s</p>\n',keywords);
end;

if is_review,
    fprintf(fout,'<p><b>This is a review with %s references.</b></p>\n',references);
end;

fprintf(fout,'<p>Published in %s language.</p>\n',language);

fprintf(fout,'<p><i>downloaded from PubMed on: %s</i></p>\n',datestr(now));
fprintf(fout,'</body>\n');
fprintf(fout,'</html>\n');
fclose(fout);