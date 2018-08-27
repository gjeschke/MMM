function ensemble=assess_ensemble(inname,outname,esize,min_GA341)
% Sorts a list of Modeller models assessed by both DOPE and GA341 score and
% provides a new list with at most esize models
% all models in the new list have a GA341 score of at least 0.75 or the
% minimum value specified by min_GA341
% the are sorted by ascending DOPE score
%
% inname    filename (path included) of the Modeller log file
% outname   filename for the extracted list of models, an empty string
%           suppresses output to file
% esize     number of models to be extracted
%
% ensemble  list of all extracted models, struct with fields
%           .name   file name of model
%           .molpdf Modeller objective function
%           .DOPE   DOPE score
%           .GA341  GA341 score
%
% Remarks:
%
% The Modeller Python script should contain an automodel object definition
% and update analogous to
%
% a = automodel(env,
%              alnfile  = 'PutP_Olkhova_2XQ2.ali',     # alignment filename
%              knowns   = '2XQ2A',              # codes of the templates
%              sequence = 'PutPOX',        # code of the target
%	      assess_methods=(assess.DOPE, assess.GA341)) # model assessment              
% a.starting_model= 1                 # index of the first model
% a.ending_model  = 100               # index of the last model
%                                     # (determines how many models to calculate)
% a.make()                            # do the actual homology modeling
%
% the model list then has the following columns
%
% Filename                          molpdf     DOPE score    GA341 score
%
% asize=0 indicates failure of homology modeling
% the output file is written only if models are encountered in the log file
% and at least one of them has GA341 >= 0.7
%
% G. Jeschke, 2011

if nargin<4,
    min_GA341=0.75;
end;

key='>> Summary of successfully produced models:';

asize=0;

fid=fopen(inname,'r');

molpdf=zeros(1,1000);
DOPE=zeros(1,1000);
GA341=zeros(1,1000);

poi=0;
if fid~=-1,
    isheader=1;
    while isheader,
        tline = fgetl(fid);
        if ~ischar(tline), break, end;
        if strfind(tline,key),
            tline = fgetl(fid); % headline
            tline = fgetl(fid); % dashes
            tline = fgetl(fid); % first model description
            while ~isempty(tline);
                nonsense=textscan(tline,'%s');
                fields=nonsense{1};
                if length(fields)<4,
                    disp('Wrong format of model list');
                    return
                end;
                poi=poi+1;
                names{poi}=fields{1};
                molpdf(poi)=str2double(strtrim(fields{2}));
                DOPE(poi)=str2double(strtrim(fields{3}));
                GA341(poi)=str2double(strtrim(fields{4}));
                tline = fgetl(fid); % first model description
            end;
        end;
    end;
end;

fclose(fid);

if poi<1,
    add_msg_board('No model found in log file');
    return;
end;

molpdf=molpdf(1:poi);
DOPE=DOPE(1:poi);
GA341=GA341(1:poi);

if max(GA341)<min_GA341,
    add_msg_board(sprintf('No model has GA341 >= %5.2f\n',min_GA341));
    add_msg_board(sprintf('Best GA341 is %5.2f\n',max(GA341)));
    ensemble=[];
    return
end;

[sDOPE,sind]=sort(DOPE);

esize0=poi;
poi=0;
k=0;

if ~isempty(outname),
    silent=false;
    outf=fopen(outname,'wt');

    fprintf(outf,'Reduced Modeller ensemble out of %i models\n\n',esize0);
    fprintf(outf,'#      Filename                   molpdf        DOPE score       GA341 score\n\n');
else
    silent=true;
end;

rej=0;
GA341sum=0;
while k<esize0 && poi<esize,
    k=k+1;
    ind=sind(k);
    if GA341(ind)>=min_GA341,
        poi=poi+1;
        ensemble(poi).name=names(ind);
        ensemble(poi).molpdf=molpdf(ind);
        ensemble(poi).DOPE=DOPE(ind);
        ensemble(poi).GA341=GA341(ind);
        if ~silent,
            fprintf(outf,'%3i    %s     %12.5f%16.5f%12.5f\n',poi,names{ind},molpdf(ind),DOPE(ind),GA341(ind));   
        end;
        GA341sum=GA341sum+GA341(ind);
    else
        rej=rej+1;
    end;
end;

if ~silent,
    fclose(outf);
end;

add_msg_board(sprintf('%i model(s) could be extracted with a mean GA341 score of %5.3f.\n',poi,GA341sum/poi));
add_msg_board(sprintf('%i model(s) were rejected because of too low GA341 score.\n',rej));
