function message=wr_BibTex(fname,references)
% Writes a reference list to a BibTex library file (.bib)
%
% G. Jeschke, 2009

message.error=0;
message.text='No error.';

fid=fopen(fname,'w');

if fid~=-1,
    for rec=1:length(references),
        fprintf(fid,'@');
        switch references(rec).type
            case {1,9}
                fprintf(fid,'article{');
            case 2
                fprintf(fid,'book{');
            case 3
                fprintf(fid,'inbook{');
            case 4
                fprintf(fid,'misc{'); % webpage
            case 5
                fprintf(fid,'conference{'); % conference
            case 6
                fprintf(fid,'phdthesis{'); % dissertation
            case 7
                fprintf(fid,'misc{'); % preprint
            case 8
                fprintf(fid,'misc{'); % patent
        end;
        fprintf(fid,'%s,\n',references(rec).short);
        if ~isempty(references(rec).authors),
            authors=textscan(references(rec).authors,'%s','delimiter',';');
            all_authors=' author    = {';
            if ~isempty(authors),
                if ~isempty(authors{1}),
                    for k=1:length(authors{1}),
                        author0=strtrim(char(authors{1}(k)));
                        author='';
                        for kk=1:length(author0),
                            if author0(kk)~='.' && author0(kk)~=',',
                                author=[author author0(kk)];
                            end;
                        end;                                
                        [surname,initials]=strtok(author);
                        initials=strtrim(initials);
                        surname=strtrim(surname);
                        if length(initials)==1,
                            initials=[initials '. '];
                        end;
                        if length(initials)==2,
                            initials=[initials(1) '. ' initials(2) '. '];
                        end;
                        if length(initials)>2,
                            initials=[initials ' '];
                        end;
                        if ~isempty(initials),
                            author=[initials surname];
                        else
                            author=surname;
                        end;
                        if k==1,
                            all_authors=[all_authors author];
                        else
                            all_authors=[all_authors ' and ' author];
                        end;
                    end;
                    fprintf(fid,'%s},\n',all_authors);
                end;
            end;
        else
            fprintf(fid,' author     = {anonymous},\n');
        end;
        fprintf(fid,' title     = {%s},\n',references(rec).title);
        if ~isempty(references(rec).publisher),
            fprintf(fid,' publisher = {%s},\n',references(rec).publisher);
        end;
        if ~isempty(references(rec).city),
            fprintf(fid,' location  = {%s},\n',references(rec).city);
        end;
        if ~isempty(references(rec).publisher) && references(rec).type==6,
            fprintf(fid,' school    = {%s},\n',references(rec).publisher);
        end;
        fprintf(fid,' journal   = {%s},\n',references(rec).journal);
        fprintf(fid,' year      = {%s},\n',references(rec).year);
        if references(rec).type~=8,
            fprintf(fid,' volume    = {%s},\n',references(rec).volume);
            fprintf(fid,' number    = {%s},\n',references(rec).issue);
        else
            fprintf(fid,' country   = {%s},\n',references(rec).volume);
            fprintf(fid,' number    = {%s},\n',references(rec).issue);
        end;
        poi=strfind(references(rec).pages,'-');
        if isempty(poi) || length(poi)>1 || poi<=1 || poi==length(references(rec).pages),
            fprintf(fid,' pages     = {%s},\n',references(rec).pages);
        else
            fprintf(fid,' pages     = {%s--%s},\n',references(rec).pages(1:poi-1),references(rec).pages(poi+1:end));
        end;
        if ~isempty(references(rec).book_title),
            fprintf(fid,' booktitle = {%s},\n',references(rec).book_title);
        end;
        if ~isempty(references(rec).chapter),
            fprintf(fid,' chapter   = {%s},\n',references(rec).chapter);
        end;        
        if ~isempty(references(rec).editors),
            fprintf(fid,' editor    = {%s},\n',references(rec).editors);
        end;
        if ~isempty(references(rec).URL),
            fprintf(fid,' URL       = {%s},\n',references(rec).URL);
        end;
%         if ~isempty(references(rec).DOI),
%             fprintf(fid,' note      = {doi: %s}\n',references(rec).DOI);
%         elseif references(rec).PMID>0,
%             fprintf(fid,' note      = {PMID: %i}\n',references(rec).PMID);
%         else
%             fprintf(fid,' note      = {no electronic access information}\n');
%         end;
        fprintf(fid,'}\n\n');
    end;
    fclose(fid);
else
    message.error=1;
    message.text='ERROR. File could not be opened for export.';
end;
