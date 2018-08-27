function compile_manual(mode)
% function compile_manual
%
% compiles a manual in PDF format via an intermediate RTF file from the set
% of MMM help files in HTML format
%
% the arrangement of the help files in the manual is specified in the file
% manual_template.txt in MMM's help subdirectory
% the title page is given in the file manual.html in MMM's help
% subdirectory (not directly accessible in online help)
% 
% the intermediate RTF file is stored as manual_compiled.rtf in MMM's help
% subdirectory
% all other intermediate files are stored in the tmp subdirectory
% the output file manual.pdf is stored in MMM's help subdirectory
%
% Input parameters:
%
% mode  can be 'html' or 'rtf' or 'pdf' or 'all' (case-insensitive) or 
%       empty, if it is empty or missing, the full compilation is performed 
%       with a query about manual editing of the intermediate RTF file
%       if 'all', there is no query 
%       if 'html', only the HTML file of the manual is compiled, this has
%       the file name manual_compiled.html
%       if 'rtf', only the RTF file is prepared
%       if 'pdf', a user-defined RTF file is compiled into the manual,
%       which allows to have an edited version (page numbers, small
%       formatting issues, manual insertions) that is not overwritten on
%       next compilation
%
% this function is made for developers, not users
% it requires the two licensed shareware programs htmltortf.exe and
% doc2pdf.exe, which are reasonably priced, but not for free
% these programs must be on the Matlab path and you must have a valid
% license for them
% you can also manually prepare the PDF version from the RTF version using
% OpenOffice or MS Word and a PDF converter of your choice (we recommend
% the fast freeware program CutePDF)
%
% G. Jeschke, 2009

global general

if nargin<1 || isempty(mode), mode=''; end;
mode=upper(mode);

clc
fclose('all');

pdf_out=fullfile(general.help_files,'manual.pdf');
outfile=fullfile(general.help_files,'manual_compiled.rtf');
current=pwd;

disp('--- MMM manual compiler ---');
disp(' ');

if ~strcmpi(mode,'PDF'),
    disp('Assembling the manual HTML file from help files');

    % First pass, generates the Table of contents (without page numbers)
    section=0;
    subsection=0;
    subsubsection=0;
    fid=[];

    ofile=[general.tmp_files 'contents.html'];
    ifile2=ofile;
    outf=fopen(ofile,'w');
    fprintf(outf,'<h1>Contents</h1><p>\n');
    section_tags=':';
    sections={};
    secpoi=0;
    while 1,
        [type,arg,fid]=get_template_line(fid);
        if fid==-1, break; end;
        %    disp(sprintf('Type: %i, Argument: |%s|',type,arg));
        switch type
            case 1
                section=section+1;
                fprintf(outf,'<h2>%i %s</h2><p>',section,arg);
                subsection=0;
                section_tags=[section_tags strtrim(arg) ':'];
                secpoi=secpoi+1;
                sections{secpoi}=sprintf('%i',section);
            case 2
                subsection=subsection+1;
                headline=get_headline(arg);
                fprintf(outf,'<h3>%i.%i %s</h3><p>',section,subsection,headline);
                subsubsection=0;
                section_tags=[section_tags strtrim(arg) ':'];
                secpoi=secpoi+1;
                sections{secpoi}=sprintf('%i.%i',section,subsection);
            case 3
                subsubsection=subsubsection+1;
                headline=get_headline(arg);
                fprintf(outf,'<h4>%i.%i.%i %s</h4><p>',section,subsection,subsubsection,headline);
                section_tags=[section_tags strtrim(arg) ':'];
                secpoi=secpoi+1;
                sections{secpoi}=sprintf('%i.%i.%i',section,subsection,subsubsection);
        end;
    end;
    fclose(outf);

    disp('First pass finished');

    % Second pass, generates the intermediate HTML file

    ofile=fullfile(general.help_files,'manual_compiled.html');
    if exist(ofile,'file'),
        delete(ofile);
    end;
    ifile1=[general.help_files 'manual.html'];

    wfid=fopen(ofile,'w');
    rfid=fopen(ifile1,'r');
    while 1,
        tline=fgetl(rfid);
        if ~ischar(tline); break; end;
        fprintf(wfid,'%s\n',tline);
    end;
    fclose(rfid);
    fprintf(wfid,'<p style="page-break-before: always">\n');
    rfid=fopen(ifile2,'r');
    while 1,
        tline=fgetl(rfid);
        if ~ischar(tline); break; end;
        fprintf(wfid,'%s\n',tline);
    end;
    fclose(rfid);

    section=0;
    subsection=0;
    subsubsection=0;
    fid=[];

    while 1,
        [type,arg,fid]=get_template_line(fid);
        if fid==-1, break; end;
%        disp(sprintf('Type: %i, Argument: |%s|',type,arg));
        switch type
            case 1
                section=section+1;
                fprintf(wfid,'<p style="page-break-before: always">\n');
                fprintf(wfid,'<h1>%i %s</h1>',section,arg);
                subsection=0;
            case 2
                subsection=subsection+1;
                % fprintf(wfid,'<p>\n');
                insert_html_file(wfid,arg,sprintf('<h2>%i.%i',section,subsection),'</h2>',sections,section_tags);
                fprintf(wfid,'</p><h3>  </h3><p>\n');
                subsubsection=0;
            case 3
                subsubsection=subsubsection+1;
                % fprintf(wfid,'<p>\n');
                insert_html_file(wfid,arg,sprintf('<h2>%i.%i.%i',section,subsection,subsubsection),'</h2>',sections,section_tags);
                fprintf(wfid,'</p><h3>  </h3><p>\n');
        end;
    end;
    fprintf(wfid,'</body>\n</html>\n');
    fclose(wfid);

    if strcmpi(mode,'HTML'),
        return
    end;

    disp('Second pass finished');
    disp('Starting HTML to RTF conversion');

    fullpath=which('htmltortf.exe');
    cmdpath=fullpath(1:end-length('htmltortf.exe'));
    cd(cmdpath);

    cmd='htmltortf.exe -image yes -align justify -kh yes -kf yes -ks yes -kc yes -psize a4 -ml 20 -mr 20 -mt 20 -mb 25';
    infile=fullfile(general.help_files,'manual_compiled.html');
    outfile=fullfile(general.help_files,'manual_compiled.rtf');
    fullcmd=[cmd ' ' infile ' ' outfile];
    [status, result] = dos(fullcmd, '-echo');


    disp(sprintf('RTF compilation finished with status: %i',status));
    disp(sprintf('RTF file %s was generated.',outfile));
    disp(' ');
else % PDF file from user-defined RTF
    currdir=pwd;
    cd(general.help_files);
    [fname,pname]=uigetfile('*.rtf','Run script file','manual_compiled.rtf');
    cd(currdir);
    if isequal(fname,0) || isequal(pname,0)
        disp('RTF conversion to PDF cancelled by user');
    else
        outfile=fullfile(pname,fname);
    end;
end;

if isempty(mode),
    disp('This is the time to insert page numbers in the Table of Contents');
    ButtonName = questdlg(sprintf('Open the RTF file "%s" with MS Word or OpenOffice and insert page numbers. Save and close RTF file (same name). Click OK when finished.',outfile), 'Table of Contents needs page numbers', 'OK', 'Skip', 'OK');
    if strcmp(ButtonName,'Skip'),
        disp('Continuing without page numbers.');
    end;
    disp(' ');
end;

if ~strcmpi(mode,'RTF'),
    disp('Starting RTF to PDF conversion');

    fullpath=which('doc2pdf-0_7_1.exe');
    cmdpath=fullpath(1:end-length('doc2pdf-0_7_1.exe'));
    cd(cmdpath);
    cmd=['doc2pdf-0_7_1 -b 24 -p 0 -i ' outfile ' -o ' pdf_out];
    [status, result] = dos(cmd);
    disp(sprintf('Status: %i',status));
    if status==0,
        disp(sprintf('PDF file %s was generated.',pdf_out));
    else
        disp('Error in PDF generation.');
    end;
end;
cd(current);

function [type,arg,fid]=get_template_line(fid)
% function [type,arg,fid]=get_template_line(fid)
%
% reads lines from the template file until a non-comment line is
% encountered
% this line is interpreted and the line type and argument are returned
% 
% fid   file identifier for the template file, if argument is empty the
%       file is opened, empty fid is returned, if file cannot be opened
% type  line type
%       0   verbatim (line is to be included directly into output HTML file)
%           the argument arg is the whole (trimmed) line
%       1   section, the argument is the trimmed section title
%       2   include HTML file as subsection, the argument is the file name
%       3   include HTML file as subsubsection, argument is the file name
%       -1  file end, file is already closed
%       empty, if file does not exist

global general

type=[];
arg=[];

if isempty(fid),
    fname=[general.help_files 'manual_template.txt'];
    fid=fopen(fname,'r');
    if fid==-1,
        fid=[];
        disp('Template file not found.');
        return
    end;
end;

comment=true;
while comment
    tline=fgetl(fid);
    if ~ischar(tline), fclose(fid); fid=-1; break, end
    tline=strtrim(tline);
    if ~strcmp(tline(1),'%'),
        comment=false;
    end;
end;

if fid==-1,
    return;
end;

type=0;
if tline(1)=='$',
    type=1;
    arg=strtrim(tline(2:end));
elseif tline(1)=='<',
    if tline(2)=='<',
        type=3;
        arg=strtrim(tline(3:end));
    else
        type=2;
        arg=strtrim(tline(2:end));
    end;
end;

function headline=get_headline(htmlfile)
% Gets the headline information of an MMM help file (<div
% class="pagetitle">)
%
% error messages # HTML file missing # or # headline missing # are supplied
% instead of the headline if necessary

headline='# unknown error #';
rfile=fopen(htmlfile,'r');
if rfile==-1,
    headline='# HTML file missing #';
    return
end;

headline='# headline missing #';
while 1,
    tline=fgetl(rfile);
    if ~ischar(tline), fclose(rfile); break, end
    tline=strtrim(tline);
    s=strfind(lower(tline),'<div class="pagetitle">');
    if ~isempty(s) && length(tline)>s(1)+length('<div class="pagetitle">'),
        tline=tline(s(1)+length('<div class="pagetitle">'):end);
        en=strfind(lower(tline),'</div>');
        if ~isempty(en) && en(1)>1,
            headline=tline(1:en(1)-1);
        else
            headline=tline;
        end;
        break
    end;
end;
fclose(rfile);

function insert_html_file(fid,htmlfile,hl_start,hl_end,sections,section_tags)
% Insert an HTML file into the file with identifier fid
% use the headline start string hL_start and headline end comman hl_end
%
% error messages # HTML file missing # or # headline missing # are inserted
% instead of the headline if necessary

rfile=fopen(htmlfile,'r');
if rfile==-1,
    headline='# HTML file missing #';
    fprintf(fid,'%s%s%s\n',hl_start,headline,hl_end);
    return
end;

headline='# headline missing #';
in_links=0;
insert=0;
open_link=[];
nl=0;
isheader=false;
while 1,
    tline=fgetl(rfile);
    nl=nl+1;
    if ~ischar(tline), break, end
    tline=strtrim(tline);
    tline_in=tline;
    if isheader && length(tline)>=3,
        if strcmpi(tline(1:3),'<p>'),
            if length(tline)>3,
                tline=tline(4:end);
            else
                isheader=false;
                continue;
            end;
        end;
    end;
    s=strfind(lower(tline),'<div class="pagetitle">');
    if ~isempty(s) && length(tline)>s(1)+length('<div class="pagetitle">'),
        tline=tline(s(1)+length('<div class="pagetitle">'):end);
        en=strfind(lower(tline),'</div>');
        if ~isempty(en) && en(1)>1,
            headline=tline(1:en(1)-1);
        else
            headline=tline;
        end;
        fprintf(fid,'%s %s%s\n',hl_start,headline,hl_end);
    end;
    % skip over the links
    s=strfind(lower(tline),'<pre class="links">');
    if ~isempty(s),
        in_links=1;
    end;
    if in_links,
        s=strfind(lower(tline),'</pre>');
        if ~isempty(s),
            insert=1;
        end;
    end;
    if insert,
        s=strfind(lower(tline),'<div class="subsubtitle">version');
        if ~isempty(s),
            insert=0;
            break
        elseif ~in_links,
            oldlinkpos=0;
            if ~isempty(open_link),
                pos=strfind(lower(tline),'</a>');
                if ~isempty(pos),
                    oldlinkpos=pos(1)+3;
                    tline0=tline;
                    tline=[tline0(1:pos(1)+3)  sprintf(' (<i>Section %s</i>) ',sections{open_link})];
                    if pos(1)+3<length(tline0),
                        tline=[tline tline0(pos(1)+4:end)];
                    end;
                    open_link=[];
                end;
            end;
            [tline,links,ends]=extract_links(tline,1);
            if ~isempty(links),
                tline0=tline;
                tline='';
                pos=1;
                for k=1:length(links),
                    if ends(k)==-1,
                        tline=[tline tline0(pos:end)];
                        id=tag2id(upper(strtrim(links{k})),upper(section_tags));
                        if ~isempty(id) && ~strcmpi(strtrim(links{k}),htmlfile),
                            open_link=id;
                        end;
                    else
                        if ends(k)<=oldlinkpos,
                            disp(sprintf('Error in HTML link to %s in line %i of file %s',links{k},nl,htmlfile));
                            disp('Previous link was not closed with </a>.');
                        end;
                        if ends(k)>length(tline0),
                            ends(k)=length(tline0);
                        end;
                        tline=[tline tline0(pos:ends(k))];
                        pos=ends(k);
                        id=tag2id(upper(strtrim(links{k})),upper(section_tags));
                        if ~isempty(id) && ~strcmpi(strtrim(links{k}),htmlfile),
                            tline=[tline sprintf(' (<i>Section %s</i>) ',sections{id})];
                        end;
                    end;
                end;
                if pos<length(tline0),
                    tline=[tline tline0(pos+1:end)];
                end;
            end;
            fprintf(fid,'%s\n',tline);
            isheader=false;
%             s=strfind(lower(tline_in),'<div class="subtitle">');
%             if ~isempty(s), isheader=true; end;
%             s=strfind(lower(tline_in),'<h1>');
%             if ~isempty(s), isheader=true; end;
%             s=strfind(lower(tline_in),'<h2>');
%             if ~isempty(s), isheader=true; end;
%             s=strfind(lower(tline_in),'<h3>');
%             if ~isempty(s), isheader=true; end;
%             s=strfind(lower(tline_in),'<h4>');
%             if ~isempty(s), isheader=true; end;
        else
            in_links=0;
        end;
    end;
end;
fclose(rfile);
if strcmp(headline,'# headline missing #'),
    fprintf(fid,'%s%s%s\n',hl_start,headline,hl_end);
end;
if insert,
    disp(sprintf('HTML file %s is missing version information',htmlfile));
end;

function [tline,links,ends]=extract_links(tline,nomarker)
% function link=extract_links(tline,mode)
% extracts links from a HTML file line tline
%
% nomarker  optional argument, if true (1), a marker part of the link (#
%           and everything after) is suppressed
% links     link file names or URLs, empty cellstring if the line does not
%           contain a link, the corresponding cell contains an empty string
%           if the link could not be extracted (HTML syntax error)
% ends      end positions of the links (last character of <\a>), for links
%           that are still open at the end of the line, -1 is returned
% tline     is cleared of leading and trailing blanks and all sequences of
%           several blanks are condensed to single blanks, positions in
%           ends refer to the cleaned up tline that is returned

links={};
ends=[];

if nargin<2,
    nomarker=0;
end;

linkstop='</a>';
lenstop=length(linkstop)-1;
tline=compact(tline);
pos=strfind(lower(tline),'<a href');
if ~isempty(pos),
    ends=-ones(1,length(pos));
    for k=1:length(pos),
        myend=-1;
        mylink='';
        in_link=0;
        past_link=0;
        completed=0;
        for kk=pos(k):length(tline),
            if in_link,
                if tline(kk)=='"',
                    in_link=0;
                    past_link=1;                    
                elseif ~completed,
                    mylink=[mylink tline(kk)];
                end;
            elseif tline(kk)=='"',
                in_link=1;
            end;
            if past_link && ~completed,
                if kk<=length(tline)-lenstop,
                    if strcmpi(tline(kk:kk+lenstop),linkstop),
                        ends(k)=kk+lenstop+1;
                        past_link=0;
                        completed=1;
                    end;
                end;
                if completed,
                    continue;
                end;
            end;
        end;
        if nomarker,
            [mylink,rem]=strtok(mylink,'#');
        end;
        links{k}=mylink;
    end;
end;