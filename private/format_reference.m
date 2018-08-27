function formatted=format_reference(reference,format,target)
% generate an output string 'formatted' from a 'reference' and a given
% citation 'format' for a certain output 'target', which can be
%
% target    'LaTeX': a Latex string
%           'bibitem': a Latex bibitem
%           'RTF': rich text format
%           'plain': plain text (the default)
%
% formats are defined in def_reference_formats.m in folder /definitions

target=lower(target);

% make author list
author={};
poi=0;
if ~isempty(reference.authors),
    authors=textscan(reference.authors,'%s','Delimiter',';');
    if ~isempty(authors{1}),
        for k=1:length(authors{1}),
            switch k,
                case 1
                    fstring=format.first_author;
                case length(authors{1})
                    if k~=1,
                        fstring=format.last_author;
                    end;
                otherwise
                    fstring=format.any_author;
            end;
            ipoi=strfind(fstring,'I');
            spoi=strfind(fstring,'S');
            if ipoi>spoi,
                istring=fstring(ipoi:end);
                sstring=fstring(1:ipoi-1);
            else
                istring=fstring(1:spoi-1);
                sstring=fstring(spoi:end);
            end;
            if ~isempty(strtrim(char(authors{1}(k)))),
                poi=poi+1;
                old_author=strtrim((char(authors{1}(k))));
                [surname,initials]=author_parts(old_author);
                spart='';
                for kk=1:length(sstring),
                    if strcmp(sstring(kk),'S'),
                        spart=[spart surname];
                    else
                        spart=[spart sstring(kk)];
                    end;
                end;
                ipart='';
                for ki=1:length(initials),
                    for kk=1:length(istring),
                        if strcmp(istring(kk),'I'),
                            ipart=[ipart initials(ki)];
                        else
                            ipart=[ipart istring(kk)];
                        end;
                    end;
                end;
                if isempty(initials),
                    spart=surname;
                    ipart='';
                end;
                if spoi>ipoi,                    
                    author{poi}=[ipart spart];
                else
                    author{poi}=[spart ipart];
                end;
            end;
        end;
    end;
end;
if isempty(author),
    author{1}='anonymus';
end;
max_authors=format.max_authors;
if isfield(format,'include_authors'),
    max_authors=format.include_authors;
end;
if max_authors<length(author),
    oa=author;
    clear author;
    for ka=1:format.max_authors,
        author{ka}=oa{ka};
    end;
    etal=1;
else
    etal=0;
end;
authors='';
astring=format.authors;
if length(author)==1,
    astring='F';
end;
for k=1:length(astring),
    switch astring(k)
        case 'F'
            authors=[authors author{1}];
            poia=k+1;
        case 'I'
            poie=k-1;
            if length(author)>2,
                for kk=2:length(author)-1,
                    if kk>2 && poie>=poia,
                    	for kkk=poia:poie,
                            authors=[authors astring(kkk)];
                        end;
                    end;
                    authors=[authors author{kk}];
                end;
                k=poi;
            end;
        case 'L'
            if length(authors)>1,
                authors=[authors author{length(author)}];
            end;
        otherwise
            authors=[authors astring(k)];
    end;
end;
% repair sequences of double delimiters
poi=findstr(authors,'; ; ');
while length(poi)>=1 && poi(1)>1,
    authors=[authors(1:poi(1)-1) authors(poi(1)+2:end)];
    poi=findstr(authors,'; ; ');
end;
poi=findstr(authors,', , ');
while length(poi)>=1 && poi(1)>1,
    authors=[authors(1:poi(1)-1) authors(poi(1)+2:end)];
    poi=findstr(authors,', , ');
end;
poi=findstr(authors,'. . ');
while length(poi)>=1 && poi(1)>1,
    authors=[authors(1:poi(1)-1) authors(poi(1)+2:end)];
    poi=findstr(authors,'. . ');
end;
if length(author)==2,
    poi=findstr(authors,', and');
    if ~isempty(poi) && poi(1)>1,
        authors=[authors(1:poi(1)-1) authors(poi(1)+1:end)];
    end;
end;

% determine first and last page
poi=strfind(reference.pages,'-');
if isempty(poi),
    spage=reference.pages;
    epage='';
elseif poi>1,
    spage=reference.pages(1:poi-1);
    if poi<length(reference.pages),
        epage=reference.pages(poi+1:end);
    else
        epage='';
    end;
else
    spage='???';
    epage='';
end;

% delimit journal words, this may overdo full stops, but otherwise we would
% need to have a list of all abbreviated titles
journal='';
if ~isempty(reference.journal),
    journal=fuzzy_full_stop(reference.journal,format.journal_delimiter);
%     jwords=textscan(reference.journal,'%s');
%     if ~isempty(jwords{1}),
%         for k=1:length(jwords{1}),
%             journal=[journal strtrim(char(jwords{1}(k))) format.journal_delimiter ' '];
%         end;
%     end;
end;
pattern=format.order;
fonts=format.font;
formatted='';
curr_font='p';
switch fonts(1),
    case 'i'
        curr_font='i';
        switch lower(target)
            case {'latex','bibitem'}
                formatted=[formatted '\textit{'];
            case 'rtf'
                formatted=[formatted '\i '];
        end;
    case 'b'
        curr_font='b';
        switch target
            case {'latex','bibitem'}
                formatted=[formatted '\textbf{'];
            case 'rtf'
                formatted=[formatted '\b '];
        end;
end;
for k=1:length(pattern),
    suppress=0;
    if ~isempty(formatted),
        if strcmp(formatted(end),'.'),
            suppress=1;
        end;
        if strcmp(formatted(end),','),
            suppress=1;
        end;
        if strcmp(formatted(end),';'),
            suppress=1;
        end;
        if strcmp(formatted(end),':'),
            suppress=1;
        end;
        if strcmp(formatted(end),'}'),
            if length(formatted)>1,
                if strcmp(formatted(end-1),'.'),
                    suppress=1;
                end;
                if strcmp(formatted(end-1),','),
                    suppress=1;
                end;
                if strcmp(formatted(end-1),';'),
                    suppress=1;
                end;
                if strcmp(formatted(end-1),':'),
                    suppress=1;
                end;
            end;
        end;
    end;
    switch pattern(k)
        case 'A'
            inpages=0;
            formatted=[formatted authors];
            if etal,
                switch target
                    case {'latex','bibitem'}
                        formatted=[formatted ' \textit{et al.}'];
                    case 'rtf'
                        formatted=[formatted ' \i et al.\i0 '];
                    case 'html'
                        formatted=[formatted ' <i> et al.</i>'];
                    otherwise
                        formatted=[formatted ' et al.'];
                end;
            end;
        case 'T'
            inpages=0;
            formatted=[formatted reference.title];
        case 'J'
            inpages=0;
            formatted=[formatted journal];
        case 'Y'
            inpages=0;
            formatted=[formatted reference.year];
        case 'V'
            inpages=0;
            formatted=[formatted reference.volume];
        case 'I'
            inpages=0;
            formatted=[formatted reference.issue];
        case 'P'
            formatted=[formatted spage];
            inpages=1;
        case 'E'
            inpages=0;
            if ~isempty(epage),
                formatted=[formatted epage];
            end;
        case '-'
            if ~isempty(epage) || ~inpages,
                switch target
                    case {'latex','bibitem'}
                        formatted=[formatted '--'];
                    otherwise
                        formatted=[formatted '-'];
                end;
            end;
        % all of the following should never occur twice in a row
        case '.'
            if ~suppress,
                formatted=[formatted '.'];
            end;
        case ','
            if ~suppress,
                formatted=[formatted ','];
            end;
        case ':'
            if ~suppress,
                formatted=[formatted ':'];
            end;
        case ';'
            if ~suppress,
                formatted=[formatted ';'];
            end;
        % space should be contracted
%         case ' '
%             if ~strcmp(formatted(end),' '),
%                 formatted=[formatted ' '];
%             end;
        otherwise
            formatted=[formatted pattern(k)];
    end;
    if k==length(pattern) || ~strcmp(fonts(k+1),curr_font),
        switch curr_font
            case 'i'
                curr_font='i';
                switch target
                    case {'latex','bibitem'}
                        formatted=[formatted '}'];
                    case 'rtf'
                        formatted=[formatted '\i0 '];
                    case 'html'
                        formatted=[formatted '</i>'];
                end;
            case 'b'
                curr_font='b';
                switch target
                    case {'latex','bibitem'}
                        formatted=[formatted '}'];
                    case 'rtf'
                        formatted=[formatted '\b0 '];
                    case 'html'
                        formatted=[formatted '</b>'];
                end;
        end;
    end;
    if k<length(pattern),
        if ~strcmp(curr_font,fonts(k+1)),
            curr_font=fonts(k+1);
            switch curr_font
                case 'i'
                    switch target
                        case {'latex','bibitem'}
                            formatted=[formatted '\textit{'];
                        case 'rtf'
                            formatted=[formatted '\i '];
                        case 'html'
                            formatted=[formatted '<i>'];
                    end;
                case 'b'
                    switch target
                        case {'latex','bibitem'}
                            formatted=[formatted '\textbf{'];
                        case 'rtf'
                            formatted=[formatted '\b '];
                        case 'html'
                            formatted=[formatted '<b>'];
                    end;
            end;
        end;
    end;
end;

if strcmp(target,'bibitem'),
    formatted=sprintf('\\bibitem{%s}\n%s\n',reference.short,formatted);
end;

function [surname,initials]=author_parts(author)
% Separates author string in format surname initials into surname and
% initials, taking care of some known composed surnames
%

surname='';
initials='';
if isempty(author); return; end;

noblesse{1}='de';
noblesse{2}='van';
noblesse{3}='von';
noblesse{4}='Mac';
noblesse{5}='Mc';
noblesse{6}='des';
noblesse{7}='dos';

for k=1:length(noblesse),
    len=length(noblesse{k});
    poi=strfind(upper(author),upper(noblesse{k}));
    for kk=1:length(poi),
        if length(author)>poi(kk)+len,
            if strcmp(author(poi(kk)+len),' '),
                author(poi(kk)+len)='~';
            end;
        end;
    end;
end;

[surname,initials]=strtok(author);

for k=1:length(surname),
    if strcmp(surname(k),'~'),
        surname(k)=' ';
    end;
end;

% defensive programming, should not happen in initials
for k=1:length(initials),
    if strcmp(initials(k),'~'),
        initials(k)=' ';
    end;
end;

surname=strtrim(surname);
initials=strtrim(initials);
initials0=initials;
poi=0;
for k=1:length(initials0);
    if ~strcmp(initials0(k),'.'),
        poi=poi+1;
        initials(poi)=initials0(k);
    end;
end;
initials=initials(1:poi);

function journal=fuzzy_full_stop(journal0,delimiter)
% adds full stops at appropriate places to abbreviated journal string if
% delimiter is a full stop, if not, the original title is returned
%
% based on National Library of Medicine Title Abbreviations rules as of
% November 17th 2009 (actually, rules from March 1st 2007)
% these rules agree (almost) with ISSN Centre rules at March 1st 2007
% see: http://www.nlm.nih.gov/pubs/factsheets/constructitle.html
%
% this is a 'fuzzy' routine, as the rules were changed in the past (and
% may change in the future) and we do not check against the (extremely
% long) abbreviation list maintained by ISSN

minlen=3;           % the minimum length for which we assume that a word may not be an abbreviation
maxlen=7;           % the maximum word length for which we assume that the word is an abbreviation
vowels='aeiouy)]';  % we assume that words ending on a vowel or a closing parentheis are 
                    % not abbreviations, there are exceptions, but they are rare

% remove all full stops that already exist
journal1=journal0;
poi=0;
for k=1:length(journal0);
    if ~strcmp(journal0(k),'.'),
        poi=poi+1;
        journal1(poi)=journal0(k);
    end;
end;
journal0=journal1(1:poi);

if ~strcmp(delimiter,'.'),
    journal=strtrim(journal0);
    return;
end;

jwords=textscan(journal0,'%s');
if ~isempty(jwords{1}),
    % One word titles are never abbreviated.
    if length(jwords{1})==1,
        journal=strtrim(char(jwords{1}(1)));
        return
    end;
    journal='';
    for k=1:length(jwords{1}),
        word=strtrim(char(jwords{1}(k)));
        if length(word)<minlen, %this is not foolproof, but pretty good
            word=[word '.'];
        elseif length(word)<=maxlen, %this is not foolproof either, longer abbreviated words exist
            if isempty(strfind(vowels,lower(word(end)))),
                word=[word '.'];
            end;
        end;
        journal=[journal word ' '];
    end;
end;

journal=strtrim(journal);