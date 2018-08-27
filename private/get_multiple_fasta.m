function alignment=get_multiple_fasta(inname)
%
% function alignment=get_multiple_fasta(inname),
%
% Reads multiple peptide (protein) sequences in FASTA format
%
% inname    filename INCLUDING extension
% alignment alignment array of struct (one for each sequence) with fields 
%           .db         data base type
%           .sequence   string that codes the sequence in single-letter format
%           .name       identifier
%           further fields depend on data base type (see Reamrks)
%
% Remarks:
% - assumed input format for the so-called FASTA defline is based on the 
%   "Sequence identifiers" section of
%   http://en.wikipedia.org/wiki/FASTA_format (as of 21.03.2011)
%
%   generally, the FASTA defline is ill-defined, but major databases use
%   well-defined subformats
%
%
% (c) G. Jeschke, 2011

    alignment=[];
    rfile=fopen(inname,'r');
    cseq=0;
    sequence=':'; % only placeholder, will be cut out at the end

    nl=0;
    while 1,
       tline=fgetl(rfile); nl=nl+1;
       if ~ischar(tline), break, end
       if ~isempty(tline),
           if char(tline(1))=='>',
               cseq=cseq+1; % current sequence number
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
               nonsense=textscan(tline(2:end),'%s','Delimiter','|');
               tags=nonsense{1};
               if isempty(tags),
                   db='unk';
               else
                   db=char(tags(1));
               end;
               switch db
                   case 'gi'
                       alignment(cseq).db='gi';
                       alignment(cseq).gi_number=char(tags(2));
                       alignment(cseq).type=char(tags(3));
                       alignment(cseq).accession=char(tags(4));
                       alignment(cseq).name=char(tags(5));
                   case 'pir'
                       alignment(cseq).db='pir';
                       alignment(cseq).name=char(tags(3));
                   case 'prf'
                       alignment(cseq).db='prf';
                       alignment(cseq).name=char(tags(3));
                   case 'sp'
                       alignment(cseq).db='sp';
                       alignment(cseq).accession=char(tags(2));
                       alignment(cseq).name=char(tags(3));
                   case 'pdb'
                       alignment(cseq).db='pdb';
                       alignment(cseq).name=char(tags(2));
                       if length(tags)>=3,
                            alignment(cseq).chain=char(tags(3));
                       else
                           alignment(cseq).chain=' ';
                       end;
                   case 'pat'
                       alignment(cseq).db='pat';
                       alignment(cseq).country=char(tags(2));
                       alignment(cseq).name=char(tags(3));
                   case 'bbs'
                       alignment(cseq).db='bbs';
                       alignment(cseq).name=char(tags(2));
                   case 'gnl'
                       alignment(cseq).db='gnl';
                       alignment(cseq).type=char(tags(2));
                       alignment(cseq).name=char(tags(3));
                   case 'ref'
                       alignment(cseq).db='ref';
                       alignment(cseq).accession=char(tags(2));
                       alignment(cseq).name=char(tags(3));
                   case 'tr'
                       alignment(cseq).db='tr';
                       alignment(cseq).accession=char(tags(2));
                       alignment(cseq).name=char(tags(3));
                   case 'lcl'
                       alignment(cseq).db='lcl';
                       alignment(cseq).name=char(tags(2));
                   case 'MMM'
                       alignment(cseq).db='MMM';
                       alignment(cseq).name=char(tags(2));
                   case 'Clustal'
                       alignment(cseq).db='Clustal';
                       alignment(cseq).name=char(tags(2));
                   case 'Modeller'
                       alignment(cseq).db='MMM'; % Modeller information is lost in FASTA format
                       alignment(cseq).name=char(tags(2));
                   otherwise
                       [type,name]=strtok(char(tags(1)),':');
                       if length(tags)>=4 && strcmpi(char(tags(2)),'PDBID')...
                               && strcmpi(char(tags(3)),'CHAIN')...
                               && strcmpi(char(tags(4)),'SEQUENCE'),
                           alignment(cseq).db='pdb';
                           [entry,chain]=strtok(char(tags(1)),':');
                           alignment(cseq).name=entry;
                           alignment(cseq).chain=chain(2:end);
                       elseif strcmpi(type,'PDB'),
                           alignment(cseq).db='pdb';
                           if length(name)>1,
                                alignment(cseq).name=name(2:end);
                                alignment(cseq).chain=' ';
                           else
                               alignment(cseq).name='unknown';
                               alignment(cseq).chain=' ';
                           end;
                       else
                           alignment(cseq).db='unk';
                           if length(tags)>1,
                               alignment(cseq).name=char(tags(1));
                           end;
                       end;
               end;
           elseif char(tline(1)) ~= ';', % neglect comments
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
    if cseq>0,
        alignment(cseq).sequence=sequence;
    end;
    fclose(rfile);

end