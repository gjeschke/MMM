function alignment=get_multiple_pir(inname)
%
% function alignment=get_multiple_pir(inname),
%
% Reads multiple peptide (protein) sequences in PIR format
%
% inname    filename INCLUDING extension
% alignment alignment array of struct (one for each sequence) with fields 
%           .sequence   string that codes the sequence in single-letter format
%           .
%
% (c) G. Jeschke, 2011

    alignment=[];

    rfile=fopen(inname,'r');
    if rfile==-1,
        return
    end;
    cseq=0;
    sequence=':'; % only placeholder, will be cut out at the end

    nl=0;
    while 1,
       tline=fgetl(rfile); nl=nl+1;
       if ~ischar(tline), break, end
       if ~isempty(tline),
           if char(tline(1))=='C',
               comment=strtrim(tline(3:end));
               add_msg_board('PIR alignment file comment:');
               add_msg_board(sprintf('%s',comment));
           end;
           if char(tline(1))=='>',
               cseq=cseq+1; % current sequence number
               seqtype=char(tline(2:3));
               idCode=char(tline(5:end));
               alignment(cseq).seqtype=seqtype;
               alignment(cseq).name=idCode;
               alignment(cseq).db='pir_native';
               if cseq>1,
                   k=strfind(sequence,'*');
                   if ~isempty(k), 
                       sequence=sequence(2:k-1); 
                   elseif length(sequence)>=2,
                       sequence=sequence(2:end);
                   else
                       sequence='';
                   end;
                   alignment(cseq-1).sequence=sequence;
                   sequence=':';
               end;
               description=fgetl(rfile); nl=nl+1;
               if ~ischar(tline), break, end
               nonsense=textscan(description,'%s','Delimiter',':');
               if ~isempty(nonsense) && length(nonsense{1})>=6, % this is recognized as Modeller format, see Modeller manual, Appendix B1
                   tags=nonsense{1};
                   alignment(cseq).db='Modeller';
                   alignment(cseq).type=char(tags(1));
                   alignment(cseq).pdb=char(tags(2));
                   alignment(cseq).first_res=char(tags(3));
                   alignment(cseq).first_chain=char(tags(4));
                   alignment(cseq).last_res=char(tags(5));
                   alignment(cseq).last_chain=char(tags(6));
                   if length(tags)==7, % more than six arguments are optional
                       alignment(cseq).description=char(tags(7));
                   end;
                   if length(tags)>=8,
                       alignment(cseq).protein=char(tags(7));
                       alignment(cseq).source=char(tags(8));
                   else
                       alignment(cseq).protein='unknown';
                       alignment(cseq).source='unknown';
                   end;
                   if length(tags)==9, % crystallographic information is optional
                       alignment(cseq).description=char(tags(9));
                   end;
                   if length(tags)>=10,
                       alignment(cseq).resolution=char(tags(9));
                       alignment(cseq).Rfactor=char(tags(10));
                   else
                       alignment(cseq).resolution='-1.0';
                       alignment(cseq).Rfactor='-1.0';
                   end;
                   if length(tags)>10,
                       alignment(cseq).description=char(tags(11));
                   else
                       alignment(cseq).description='Modeller';
                   end;
               else
                   alignment(cseq).description=description;
               end;
           elseif char(tline(1)) ~= 'C', % neglect comments
               sequence=strcat(sequence,char(tline)); % append line to sequence
           end;
       end;
    end;
    k=strfind(sequence,'*');
    if ~isempty(k), 
       sequence=sequence(2:k-1); 
    elseif length(sequence)>=2,
       sequence=sequence(2:end);
    else
       sequence='';
    end;
    alignment(cseq).sequence=sequence;
    fclose(rfile);
    
end