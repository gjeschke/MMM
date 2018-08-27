function references=rd_SciFinder(fname,global_num)
% Reads a reference list from a CAS SciFinder-exported file
% the file must have been exported in the tagged format (*.txt)
%
% Journal Articles (type 1), Books (2), Book Sections (3), and Web
% pages (4), Conference (5) and Dissertation (6) and Preprint (7) and Patent (8) are 
% explicit types, all other source types have type 0, but
% fields are still imported for them, note that CAS do not refer to webpages
% and that there is no support of digital object identifiers (DOI) and no
% real support of URLs
% note also that SciFinder does not provide Publisher Information, hence
% book references are incomplete
% the short name is set to a running global reference number 
% in square brackets, which is initialized by global_num
%
% G. Jeschke, 2009

references=[];

fid=fopen(fname,'r');
nl=0;
rec=0;
in_record=0;
school='';
inventor='';
if fid~=-1,
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        nl=nl+1;
        poi=strfind(upper(tline),'START_RECORD');
        if ~isempty(poi),
            if rec>0 && references(rec).type==6 && ~isempty(school),
                references(rec).publisher=school;
            end;
            school='';
            if rec>0 && references(rec).type==8 && ~isempty(inventor),
                references(rec).authors=inventor;
            end;
            inventor='';
            rec=rec+1;
            global_num=global_num+1;
            references(rec).format=2;
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
            in_record=1;
            continue;
        end;
        poi=strfind(upper(tline),'END_RECORD');
        if ~isempty(poi),
            in_record=0;
        end;
        if ~in_record,
            continue;
        end;
        poi=strfind(upper(tline),'FIELD');
        if isempty(poi),
            continue;
        else
            tline=tline(6:end);
        end;
        [token,arg]=strtok(tline,':');
        switch strtrim(token)
            case 'Document Type'
                type=strtrim(arg(2:end));
                if strfind(type,'Journal'),
                    references(rec).type=1;
                end;
                if strfind(type,'Book'),
                    references(rec).type=2;
                end;
                if strfind(type,'Book Section'),
                    references(rec).type=3;
                end;
                if strfind(type,'Conference'),
                    references(rec).type=5;
                end;
                if strfind(type,'Dissertation'),
                    references(rec).type=6;
                end;
                if strfind(type,'Preprint'),
                    references(rec).type=7;
                end;
                if strfind(type,'Patent'),
                    references(rec).type=8;
                end;
            case 'Journal Title',
                references(rec).journal=strtrim(arg(2:end));
            case 'Title'
                references(rec).title=strtrim(arg(2:end));
            case 'Editor'
                references(rec).editors=strtrim(arg(2:end));
            case 'Author'
                references(rec).authors=process_authors(strtrim(arg(2:end)));
            case 'Chapter'
                references(rec).chapter=strtrim(arg(2:end));
            case 'Publication Year'
                references(rec).year=strtrim(arg(2:end));
            case 'Volume'
                references(rec).volume=strtrim(arg(2:end));
            case 'Patent Country'
                if references(rec).type==8,
                    references(rec).volume=strtrim(arg(2:end));
                end;
            case 'Issue'
                references(rec).issue=strtrim(arg(2:end));
            case 'Patent Number'
                if references(rec).type==8,
                    references(rec).issue=strtrim(arg(2:end));
                end;
            case 'Page'
                references(rec).pages=strtrim(arg(2:end));
            case 'DOI'
                references(rec).DOI=strtrim(arg(2:end));
            case 'URL'
                references(rec).URL=strtrim(arg(2:end));
            case 'City'
                references(rec).city=strtrim(arg(2:end));
            case 'Corporate Source'
                school=strtrim(arg(2:end));
            case 'Inventor Name'
                inventor=strtrim(arg(2:end));
        end;
    end
    % the following would be nice, but requires a computationally
    % inefficient test of whether this short name is already used
%     if ~isempty(references(rec).year) && ~isempty(references(rec).authors),
%         first_author=strtok(references(rec).authors);
%         references(rec).short=sprintf('%s:%i',first_author,references(rec).year);
%     end;
    if rec>0 && references(rec).type==6 && ~isempty(school),
        references(rec).publisher=school;
    end;
    if rec>0 && references(rec).type==8 && ~isempty(inventor),
        references(rec).authors=inventor;
    end;
    fclose(fid);
end;

function cleaned=process_authors(author_names)
% cleans up an author list by initialing first names and removing commas 

cleaned='';
nonsense=textscan(strtrim(author_names),'%s','Delimiter',';');
authors=nonsense{1};
for k=1:length(authors),
    author=char(authors(k));
    for kk=1:length(author),
        if author(kk)==',' || author(kk)=='.',
            author(kk)=' ';
        end;
    end;
    [lastname,firstnames]=strtok(author);
    cleaned=[cleaned lastname ' '];
    firstnames=strtrim(firstnames(2:end));
    if ~isempty(firstnames),
        nonsense=textscan(firstnames,'%s','Delimiter',' ');
        first_name_list=nonsense{1};
        for kk=1:length(first_name_list),
            name=char(first_name_list(kk));
            if ~isempty(name),
                cleaned=[cleaned name(1)];
            end;
        end;
        cleaned=[cleaned ';'];
    end;
end;
if length(cleaned)>1,
    cleaned=cleaned(1:end-1);
end;