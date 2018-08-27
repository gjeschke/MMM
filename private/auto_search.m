function auto_search

% Check whether defined and active auto searches are up to date
% if not, update these searches

global model
global queries

if isfield(model,'autosearch'),
    if ~isempty(model.autosearch),
        for k=1:length(model.autosearch),
            search=model.autosearch(k);
            if search.active,
                last_update=search.date;
                days=ceil(now-last_update);
                if days>queries.update,
                    days=200;
                    add_msg_board(sprintf('Reference search %s:',search.name));
                    status=do_search(search,days+1);
                    if status,
                        search.date=now;
                        model.autosearch(k)=search;
                        update_references;
                    end;
                end;
            end;
        end;
    end;
end;

function status=do_search(search,days)
% the search terms and modes are defined in struct search
% the time since the last search is given in days

global web_adr
global queries
global general

search_string='?term=';

term_mode='[all]';
switch search.mode
    case 2
        term_mode='[tiab]';
    case 3
        term_mode='[mh]';
    case 4
        term_mode='[majr]';
end;
terms=compact(search.search_terms);
poi=strfind(terms,';');
if ~isempty(poi),
    new_terms='';
    poi0=1;
    for k=1:length(poi),
        new_terms=[new_terms terms(poi0:poi(k)-1) term_mode '+AND+'];
        poi0=poi(k)+1;
    end;
    terms=[new_terms terms(poi0:end)];
end;
if ~isempty(terms),
    search_string=[search_string,terms,term_mode];
end;

author_mode='[au]';
switch search.author_mode
    case 2
        author_mode='[1au]';
    case 3
        author_mode='[lastau]';
end;
authors=search.authors;
for k=1:length(authors),
    if strcmp(authors(k),'.'),
        authors(k)=' ';
    end;
    if strcmp(authors(k),','),
        authors(k)=' ';
    end;
end;
poi=strfind(authors,';');
if ~isempty(poi),
    new_authors='';
    poi0=1;
    for k=1:length(poi),
        new_authors=[new_authors authors(poi0:poi(k)-1) '[au]+AND+'];
        poi0=poi(k)+1;
    end;
    authors=[new_authors authors(poi0:end)];
end;
authors=compact(authors);
if ~isempty(authors),
    if length(search_string)>length('?terms='),
        search_string=[search_string '+AND+'];
    end;
    search_string=[search_string,authors,author_mode];
end;

journal=search.journal;
for k=1:length(journal),
    if strcmp(journal(k),'.'),
        journal(k)=' ';
    end;
end;
journal=compact(journal);
if ~isempty(journal),
    if length(search_string)>length('?terms='),
        search_string=[search_string '+AND+'];
    end;
    search_string=[search_string,journal,'[ta]'];
end;

search_string=[search_string,sprintf('"last+%i+days"[edat]',days)];

if search.reviews,
    if length(search_string)>length('?terms='),
        search_string=[search_string '+AND+'];
    end;
    search_string=[search_string,'Review[pt]'];
end;
switch search.records % this is bad programming, as information from the GUI is duplicated here
                        % but it does not matter much for this purpose
    case 1
        dispmax='&dispmax=10';
    case 2
        dispmax='&dispmax=20';
    case 3
        dispmax='&dispmax=50';
    case 4
        dispmax='&dispmax=100';
end;

for k=1:length(search_string),
    if strcmp(search_string(k),' '),
        search_string(k)='+';
    end;
end;

query=sprintf('%s%s%s%s',web_adr.PubMed,search_string,queries.search_PubMed,dispmax,queries.tool,queries.email);
fname=strcat(general.tmp_files,queries.PubMed_file);
[f,status]=urlwrite(query,fname);
if status==0,
    add_msg_board('ERROR: PubMed access failed in automatic search.');
    add_msg_board('Please check internet connection.');
    add_msg_board('References will be updated when model is opened next time.');
end;

function newstr=unblank(str)
% remove all blanks from an input string

newstr='';
for k=1:length(str),
    if ~strcmp(str(k),' '),
        newstr=[newstr str(k)];
    end;
end;

function newstr=compact(str)
% condense all sequences of blanks in an input string to a single blank
% and remove all leading and trailing blanks

newstr='';
str=strtrim(str);
blanks=0;
for k=1:length(str),
    if strcmp(str(k),' '),
        blanks=blanks+1;
    else
        blanks=0;
    end;
    if blanks<2,
        newstr=[newstr str(k)];
    end;
end;

function update_references

global general
global queries
global model

if ~isfield(model,'references'),
    model.references=[];
end;
fname=strcat(general.tmp_files,queries.PubMed_file);
if isempty(model.references),
    poi=0;
else
    poi=length(model.references);
end;
references=rd_MEDLINE(fname,poi);
if isempty(references),
    add_msg_board('PubMed search does not return any references.');
else
    m=length(references);
    add_msg_board(sprintf('PubMed search returned %i reference(s).',m));
    if m>0,
        if isempty(model.references),
            model.references=references;
        else
            poi=length(model.references);
            newrefs=0;
            oldrefs=0;
            for k=1:m, % add only new references, as judged by PubMed ID (if available)
                found=0;
                if references(k).PMID~=0,                   
                    for kk=1:poi,
                        if model.references(kk).PMID==references(k).PMID,
                            found=1; break;
                        end;
                    end;
                end;
                if ~found,
                    poi=poi+1;
                    newrefs=newrefs+1;
                    model.references(poi)=references(k);
                else
                    oldrefs=oldrefs+1;
                    references(k).short=model.references(kk).short;
                    if isempty(references(k).DOI),
                        references(k).DOI=model.references(kk).DOI;
                    end;
                    model.references(kk)=references(k);
                end;
            end;
        end;
        add_msg_board(sprintf('Imported %i reference(s) with new or without any PubMed IDs',newrefs)); 
        add_msg_board(sprintf('Updated %i reference(s) with already existing PubMed IDs',oldrefs)); 
    end;
end;


   
