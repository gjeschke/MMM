function references=rd_EndNote(fname,global_num)
% Reads a reference list from an EndNote-exported file
% the file must have been exported as File type: Text file (*.txt) with
% Output style: Show All Fields
% tested with EndNote X2
%
% only Journal Articles (type 1), Books (2), Book Sections (3), and Web
% pages (4) are explicitly types, all other source types have type 0, but
% fields are still imported for them
% the short name is set to a running global reference number 
% in square brackets, which is initialized by global_num
%
% G. Jeschke, 2009

references=[];

fid=fopen(fname,'r');
nl=0;
rec=0;
if fid~=-1,
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        nl=nl+1;
        if nl==1,
            tline=tline(4:end);
        end;
        [token,arg]=strtok(tline,':');
        switch token
            case 'Reference Type'
                rec=rec+1;
                global_num=global_num+1;
                references(rec).format=3;
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
                type=strtrim(arg(2:end));
                switch type
                    case 'Journal Article'
                        references(rec).type=1;
                    case 'Book'
                        references(rec).type=2;
                    case 'Book Section'
                        references(rec).type=3;
                    case 'Web Page'
                        references(rec).type=4;
                end;
            case 'Title'
                references(rec).title=strtrim(arg(2:end));
            case 'Book Title'
                references(rec).book_title=strtrim(arg(2:end));
            case 'Editor'
                references(rec).editors=strtrim(arg(2:end));
            case 'Publisher'
                references(rec).publisher=strtrim(arg(2:end));
            case 'Author'
                references(rec).authors=strtrim(arg(2:end));
            case 'Chapter'
                references(rec).chapter=strtrim(arg(2:end));
            case 'Journal'
                references(rec).journal=strtrim(arg(2:end));
            case 'Year'
                references(rec).year=strtrim(arg(2:end));
            case 'Volume'
                references(rec).volume=strtrim(arg(2:end));
            case 'Issue'
                references(rec).issue=strtrim(arg(2:end));
            case 'Pages'
                references(rec).pages=strtrim(arg(2:end));
            case 'DOI'
                references(rec).DOI=strtrim(arg(2:end));
            case 'URL'
                references(rec).URL=strtrim(arg(2:end));
            case 'City'
                references(rec).city=strtrim(arg(2:end));
        end;
    end
    fclose(fid);
end;
