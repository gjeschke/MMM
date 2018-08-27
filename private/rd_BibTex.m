function references=rd_BibTex(fname,global_num)
% Reads a reference list from an BibTex file
%
% only Journal Articles (@article, type 1), Books (@book, 2), Book Sections
% (@inbook, 3), and Web pages (@misc, 4) are explicitly 
% typed, all other source types have type 0, but
% fields are still imported for them
% the short name is initialized with the tag in the type line, the
% recognized format is as follows:
%
% @book{Abragam:89,
%  author    = {Anatole Abragam},
%  title     = {{Principles of Nuclear Magnetism}},
%  publisher = {Clarendon Press},
%  address	 = {Oxford},
%  year      = 1989
% }
%
% with the label in the first line of the record and delimited by a comma 
% and each other field in one line, other formats, although they may work 
% in LaTex are not recognized
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
        opener= strfind(tline,'@');
        if ~isempty(opener),
            rec=rec+1;
            global_num=global_num+1;
            references(rec).format=4;
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
            [type,rem]=strtok(tline(opener+1:end),'{');
            switch type
                case 'article'
                    references(rec).type=1;
                case 'book'
                    references(rec).type=2;
                case 'inbook'
                    references(rec).type=3;
                case 'misc'
                    references(rec).type=4;
            end;
            if length(rem)>=3,
                references(rec).short=strtok(rem(2:end),','); 
            else
                references(rec).short='?';
            end;
        else
            [token,arg0]=strtok(tline,'=');
            arg='';
            for k=2:length(arg0),
                if arg0(k)~='{' && arg0(k)~='}',
                    arg=[arg arg0(k)];
                end;
            end;
            for k=1:length(arg),
                if strcmp(arg(k),'~'),
                    arg(k)=' ';
                end;
            end;
            arg=strtrim(arg);
            token=strtrim(token);
            switch token
                case 'title'
                    if references(rec).type==3,
                        references(rec).book_title=no_comma(arg);
                    else
                        references(rec).title=no_comma(arg);
                    end;
                case 'chapter'
                    if references(rec).type==3,
                        references(rec).title=no_comma(arg);
                    end;
                case 'editor'
                    references(rec).editors=no_comma(arg);
                case 'publisher'
                    references(rec).publisher=no_comma(arg);
                case 'journal'
                    references(rec).journal=no_comma(arg);
                case 'author'
                    authors0=arg;
                    sep=findstr(authors0,'and');
                    poi=1;
                    authors={};
                    if ~isempty(sep),
                        for k=1:length(sep),
                            if sep(k)>poi;
                                authors{k}=strtrim([authors0(poi:sep(k)-1)]);
                                poi=sep(k)+3;
                            end;
                        end;
                        if poi<length(authors0),
                            authors{length(authors)+1}=authors0(poi:end);
                        end;
                    else
                        authors{1}=authors0;
                    end;
                    new_authors='';
                    for k=1:length(authors),
                        sep=findstr(authors{k},',');
                        if ~isempty(sep) && sep(1)<length(authors{k}),
                            initials=strtrim(authors{k}(sep(1)+1:end));
                            if sep(1)>1,
                                surname=strtrim(authors{k}(1:sep(1)-1));
                            else
                                surname='';
                            end;
                        else
                            intermed=textscan(authors{k},'%s');
                            if ~isempty(intermed),
                                intermed=intermed{1};
                                if ~isempty(intermed),
                                    surname=char(intermed(length(intermed)));
                                    initials='';
                                    if length(intermed)>1,
                                        for k=1:length(intermed)-1,
                                            initials=[initials ' ' char(intermed(k))];
                                        end
                                    end;
                                end;
                            else
                                initials='';
                                surname=authors{k};
                            end;
                        end;
                        for kk=1:length(initials),
                            if strcmp(initials(kk),'.'),
                                initials(kk)=' ';
                            end;
                        end;
                        initials=unblank(initials);
                        % keyboard
                        authors{k}=[surname ' ' initials];
                        if k==1 || isempty(new_authors),
                            new_authors=authors{k};
                        else
                            new_authors=[new_authors '; ' authors{k}];
                        end;
                    end;
                    references(rec).authors=no_comma(new_authors);
                case 'year'
                    references(rec).year=strtrim(no_comma(arg));
                case 'volume'
                    references(rec).volume=strtrim(no_comma(arg));
                case 'month'
                    references(rec).issue=strtrim(no_comma(arg));
                case 'pages'
                    references(rec).pages=no_double_hyphen(no_comma(arg));
                case 'doi'
                    references(rec).DOI=no_comma(arg);
                case 'url'
                    references(rec).URL=no_comma(arg);
                case 'address'
                    references(rec).city=no_comma(arg);
            end;
        end;
    end
    fclose(fid);
end;

function str=no_comma(str0)

str='';
for k=1:length(str0),
    if ~strcmp(str0(k),','),
        str=[str str0(k)];
    end;
end;

function str=no_double_hyphen(str0)

str='';
hyph=0;
for k=1:length(str0),
    if strcmp(str0(k),'-'),
        if ~hyph,
            hyph=1;
            str=[str str0(k)];
        end;
    else
        str=[str str0(k)];
        hyph=0;
    end;
end;