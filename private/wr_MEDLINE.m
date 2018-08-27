function message=wr_MEDLINE(fname,references)
% Writes a reference list to a MEDLINE-formatted file
%
%
% G. Jeschke, 2009

message.error=0;
message.text='No error.';

fid=fopen(fname,'w');

if fid~=-1,
    for rec=1:length(references),
        fprintf(fid,'PMID- %i\n',references(rec).PMID);
        fprintf(fid,'DP  - %s\n',references(rec).year);
        fprintf(fid,'VI  - %s\n',references(rec).volume);
        fprintf(fid,'IP  - %s\n',references(rec).issue);
        fprintf(fid,'PG  - %s\n',references(rec).pages);
        fprintf(fid,'TI  - %s\n',references(rec).title);
        if ~isempty(references(rec).authors),
            authors=textscan(references(rec).authors,'%s','delimiter',';');
            if ~isempty(authors),
                if ~isempty(authors{1}),
                    for k=1:length(authors{1}),
                        author=strtrim(char(authors{1}(k)));
                        fprintf(fid,'AU  - %s\n',author);
                    end;
                end;
            end;
        else
            fprintf(fid,'AU  - anonymous\n');
        end;
        switch references(rec).type
            case 1
                fprintf(fid,'PT  - Journal Article\n');
            case 2
                fprintf(fid,'PT  - Book\n');
            case 3
                fprintf(fid,'PT  - Book section\n');
            case 4
                fprintf(fid,'PT  - Web page\n');
            case 5
                fprintf(fid,'PT  - Conference\n');
            case 6
                fprintf(fid,'PT  - Dissertation\n');
            case 8
                fprintf(fid,'PT  - Patent\n');
            case 9
                fprintf(fid,'PT  - Review\n');
        end;
        fprintf(fid,'TA  - %s\n',references(rec).journal);
        if ~isempty(references(rec).book_title),
            fprintf(fid,'GN  - in: %s\n',references(rec).book_title);
        end;
        if ~isempty(references(rec).chapter),
            fprintf(fid,'GN  - chapter: %s\n',references(rec).chapter);
        end;        
        if ~isempty(references(rec).editors),
            fprintf(fid,'GN  - edited by: %s\n',references(rec).editors);
        end;
        if ~isempty(references(rec).publisher),
            fprintf(fid,'GN  - publisher: %s\n',references(rec).publisher);
        end;
        fprintf(fid,'EDAT- %s\n',references(rec).year);
        if ~isempty(references(rec).DOI),
            fprintf(fid,'AID - %s [doi]\n',references(rec).DOI);
        end;
        fprintf(fid,'\n');
    end;
    fclose(fid);
else
    message.error=1;
    message.text='ERROR. File could not be opened for export.';
end;
