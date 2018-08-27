function sequence=get_pir(inname,seqnum)
%
% function sequence=get_pir(inname,seqnum),
%
% Reads a peptide (protein) sequence in FASTA format
%
% inname   is the filename INCLUDING extension
% seqnum   is an optional parameter that gives the number of the sequence in
%          the PIR file, it defaults to 1
% sequence is a string or cell string that codes the sequence in single-letter format
%
% (c) G. Jeschke, 2011

    if exist('onCleanup','class')==8, % this is done only for Matlab version 
                                      % 2008a and later, where onCleanup exists
        c = onCleanup(@myCleanup);
    end;

    if nargin<2, seqnum=1; end;

    rfile=fopen(inname,'r');
    cseq=0;
    sequence=':'; % only placeholder, will be cut out at the end

    nl=0;
    while 1
       tline = fgetl(rfile);
       if ~ischar(tline), break, end
       nl=nl+1;
       if ~isempty(tline),
           if char(tline(1))=='>',
               cseq=cseq+1; % current sequence number
               tline=fgetl(rfile); nl=nl+1;
           elseif char(tline(1)) ~= 'C' && cseq==seqnum, % neglect comments
               k=strfind(tline,'*');
               if ~isempty(k), tline=tline(1:k-1); end;
               sequence=strcat(sequence,tline); % append line to sequence
           end;
       end;
       if cseq>seqnum, break; end; % if we are already past the wanted sequence, stop
    end;
    sequence=sequence(2:length(sequence));
    fclose(rfile);
    
    function myCleanup

        if ~isempty(rfile==fopen('all')),
            fclose(rfile);
        end;
    end

end