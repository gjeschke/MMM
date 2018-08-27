function message=wr_ISI_CGI(fname,references)
% Writes a reference list to a ISI CGI-formatted file for EndNote import
% use Import option ISI-CE in EndNote
% tested with EndNote X2
%
% G. Jeschke, 2009

message.error=0;
message.text='No error.';

fid=fopen(fname,'w');

if fid~=-1,
    fprintf(fid,'FN ISI Export Format\nVR 1.0\n');
    for rec=1:length(references),
        switch references(rec).type
            case {1,9}
                fprintf(fid,'PT J\n');
            case 2
                fprintf(fid,'PT B\n');
            case 3
                fprintf(fid,'PT S\n');
            case 4
                fprintf(fid,'PT W\n');
            case 5
                fprintf(fid,'PT C\n');
            case 6
                fprintf(fid,'PT D\n');
            case 8
                fprintf(fid,'PT P\n');
        end;
        if ~isempty(references(rec).authors),
            authors=textscan(references(rec).authors,'%s','delimiter',';');
            if ~isempty(authors),
                if ~isempty(authors{1}),
                    for k=1:length(authors{1}),
                        author=strtrim(char(authors{1}(k)));
                        if k==1,
                            author=['AU ' author];
                        else
                            author=['   ' author];
                        end;
                        fprintf(fid,'%s\n',author);
                    end;
                end;
            end;
        else
            fprintf(fid,'AU anonymous\n');
        end;
        fprintf(fid,'TI %s\n',references(rec).title);
        if references(rec).type==1,
                fprintf(fid,'DT Article\n');
        end;
        if ~isempty(references(rec).publisher),
            fprintf(fid,'PU %s\n',references(rec).publisher);
        end;
        if ~isempty(references(rec).city),
            fprintf(fid,'PI %s\n',references(rec).city);
        end;
        fprintf(fid,'JI %s\n',references(rec).journal);
        fprintf(fid,'PY %s\n',references(rec).year);
        fprintf(fid,'VL %s\n',references(rec).volume);
        if references(rec).type~=8,
            fprintf(fid,'IS %s\n',references(rec).issue);
        else
            fprintf(fid,'PN %s\n',[references(rec).volume references(rec).issue]);
        end;
        if references(rec).type==1,
            fprintf(fid,'DT Article\n');
        elseif references(rec).type==9,
            fprintf(fid,'DT Review\n');
        end;
        poi=strfind(references(rec).pages,'-');
        if isempty(poi) || length(poi)>1 || poi<=1 || poi==length(references),
            fprintf(fid,'BP %s\n',references(rec).pages);
        else
            fprintf(fid,'BP %s\n',references(rec).pages(1:poi-1));
            fprintf(fid,'EP %s\n',references(rec).pages(poi+1:end));
        end;
        if ~isempty(references(rec).book_title),
            fprintf(fid,'SE %s\n',references(rec).book_title);
        end;
        if ~isempty(references(rec).chapter),
            fprintf(fid,'BS %s\n',['Chapter ' references(rec).chapter]);
        end;        
        if ~isempty(references(rec).editors),
            fprintf(fid,'ED %s\n',references(rec).editors);
        end;
        if ~isempty(references(rec).URL),
            fprintf(fid,'UR %s\n',references(rec).URL);
        end;
        if ~isempty(references(rec).DOI),
            fprintf(fid,'DI %s\n',references(rec).DOI);
        end;
        fprintf(fid,'ER\n');
        fprintf(fid,'  \n');
    end;
    fprintf(fid,'EF\n');
    fclose(fid);
else
    message.error=1;
    message.text='ERROR. File could not be opened for export.';
end;
