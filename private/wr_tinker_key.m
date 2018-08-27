function msg = wr_tinker_key(basname,options,selection)
% msg = wr_tinker_key(options,selection)
%
% Writes a Tinker keyword control file as per Section 7 of the Tinker User
% Guide, only part of Tinker functionality is supported
%
% basname   name of the input molecular system file, if empty, a default
%           file tinker.key is written
% options   option declaration, see source code for functionality and
%           defaults
% selection list of atom numbers selected as active (default) or inactive
%
% msg       message on success, msg.error is error code (0 for success) and
%           msg.text the error message
%
% G. Jeschke, 12.7.2017

msg.error = 0;
msg.text = 'Tinker keyword control file written';

if ~exist('basname','var') || isempty(basname)
    basname = 'tinker';
end

if ~exist('options','var') || isempty(options)
    options.active = true;
end

if ~exist('selection','var') || isempty(selection)
    selection = [];
end

% input/output files will be located one directory above the parameters directory

tinker_path = get_tinker_path;
if isempty(tinker_path)
    msg.error = 1;
    msg.text = 'Tinker distribution (file amber99.prm) not found on Matlab path.';
end

keyfile = fullfile(tinker_path,strcat(basname,'.key'));
fid = fopen(keyfile,'wt');
if fid == -1
    msg.error = 2;
    msg.text = 'Tinker keyword control file could not be written.';
    return
end

if isfield(options,'forcefield')
    fprintf(fid,'parameters %s\n',fullfile('params',options.forcefield));
else
    fprintf(fid,'parameters %s\n',fullfile('params','amber99'));
end

if ~isempty(selection)
    if options.active
        fprintf(fid,'active ');
    else
        fprintf(fid,'inactive ');
    end
    selection = sort(selection);
    selection = selection(selection > 0); % excludes any selected atoms that do not exist in the molecule
    block = false;
    poi = 1;
    while poi < length(selection)
        if block
            if selection(poi+1) ~= selection(poi)+1
                block = false;
                fprintf(fid,' %i',selection(poi));
            end
        else
            if selection(poi+1) == selection(poi)+1
                block = true;
                fprintf(fid,' -%i',selection(poi));
            else
                fprintf(fid,' %i',selection(poi));
            end
        end
        poi = poi + 1;
    end
    fprintf(fid,' %i\n',selection(poi));
end

if isfield(options,'solvation')
    fprintf(fid,'solvate %s\n',options.solvation);
end

if isfield(options,'restrain_d')
    [nr,m] = size(options.restrain_d);
    if m < 2
        msg.error = 3;
        msg.text = 'Distance restraint matrix has less than two columns.';
        fclose(fid);
        return
    end
    for k = 1:nr
        fprintf(fid,'restrain-distance %i %i',options.restrain_d(k,1:2));
        switch m
            case 3
                fprintf(fid,' %5.1f',options.restrain_d(k,3));
            case 5
                fprintf(fid,' %5.1f %5.2f %5.2f',options.restrain_d(k,3:5));
            otherwise
                msg.error = 4;
                msg.text = 'Distance restraint matrix must have 2, 3 or 5 columns.';
                fclose(fid);
                return         
        end
        fprintf(fid,'\n');
    end
end


if isfield(options,'restrain_t')
    [nr,m] = size(options.restrain_t);
    if m < 4
        msg.error = 5;
        msg.text = 'Torsion restraint matrix has less than four columns.';
        fclose(fid);
        return
    end
    for k = 1:nr
        fprintf(fid,'restrain-torsion %i %i %i %i',options.restrain_t(k,1:4));
        switch m
            case 5
                fprintf(fid,' %5.1f',options.restrain_t(k,5));
            case 7
                fprintf(fid,' %5.1f %5.2f %5.2f',options.restrain_t(k,5:7));
            otherwise
                msg.error = 6;
                msg.text = 'Torsion restraint matrix must have 4, 5 or 7 columns.';
                fclose(fid);
                return         
        end
        fprintf(fid,'\n');
    end
end

if isfield(options,'restrain_a')
    [nr,m] = size(options.restrain_a);
    if m < 4
        msg.error = 7;
        msg.text = 'Angle restraint matrix has less than three columns.';
        fclose(fid);
        return
    end
    for k = 1:nr
        fprintf(fid,'restrain-angle %i %i %i',options.restrain_a(k,1:3));
        switch m
            case 4
                fprintf(fid,' %5.1f',options.restrain_a(k,4));
            case 6
                fprintf(fid,' %5.1f %5.2f %5.2f',options.restrain_a(k,4:6));
            otherwise
                msg.error = 8;
                msg.text = 'Angle restraint matrix must have 3, 4 or 6 columns.';
                fclose(fid);
                return         
        end
        fprintf(fid,'\n');
    end
end

if isfield(options,'stepmax')
    fprintf(fid,'stepmax %5.2f\n',options.stepmax);
end

fclose(fid);


