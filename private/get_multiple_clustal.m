function alignment=get_multiple_clustal(inname)
%
% function alignment=get_multiple_clustal(inname),
%
% Reads multiple peptide (protein) sequences in CLUSTAL W format
%
% inname    filename INCLUDING extension
% alignment alignment array of struct (one for each sequence) with fields 
%           .sequence   string that codes the sequence in single-letter format
%           .name       sequence name
%
% (c) G. Jeschke, 2011


    alignment=[];
    
    rfile=fopen(inname,'r');

    taglist='\';
    nl=0;
    while 1,
       tline=fgetl(rfile); nl=nl+1;
       if ~ischar(tline), break, end
       if nl==1,
           if ~strcmpi(tline(1:7),'CLUSTAL') && ~strcmpi(tline(1:6),'MUSCLE'),
               break 
           end;
           tline=fgetl(rfile); nl=nl+1;
           if ~ischar(tline), break, end
       end;
       if ~isempty(strtrim(tline)),
           nonsense=textscan(tline,'%s');
           if char(tline(1))~=' ' && ~isempty(nonsense) && length(nonsense{1})>=2, % line with sequence information
               tags=nonsense{1};
               id=tag2id(char(tags(1)),taglist,[],'\'); % find sequence number
               if isempty(id),
                   taglist=[taglist char(tags(1)) '\']; % add tag, if sequence was still unknown
                   id=tag2id(char(tags(1)),taglist,[],'\');
                   alignment(id).db='Clustal';
                   alignment(id).name=char(tags(1));
                   alignment(id).sequence='';
               end;
               alignment(id).sequence=[alignment(id).sequence char(tags(2))];
           end;
       end;
    end;
    fclose(rfile);
    

end