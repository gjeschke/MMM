function [net_grad,gradients,msg] = get_Tinker_gradient(full_bas,full_force_field,correspondence)

msg.error=1;
msg.text='No result line';

net_grad=[];
gradients=[];

cmd=sprintf('testgrad %s %s Y N N',full_bas,full_force_field);

dospath=which('testgrad.exe');
if isempty(dospath),
    msg.error=2;
    msg.text='Tinker software not found on Matlab path.';
    add_msg_board('This feature requires Tinker''s testgrad from the Ponder lab');
    add_msg_board('ERROR: Tinker could not be found on the Matlab path.');
    add_msg_board('Please check whether Tinker is installed and the path set.');
    % add_msg_board('(see also help browser)');
    % webcall(entry,'-helpbrowser');
    return
end;
[modpath, modcmd] = fileparts(dospath);
my_dir=pwd;
cd(modpath);
tic;
[s, w] = dos(cmd);
runtime=toc;
add_msg_board(sprintf('Tinker testgrad was running %5.1f s',runtime));
if s~=0,
    rem=w;
    while ~isempty(rem),
        [token,rem]=strtok(rem,char(10));
        if ~isempty(token),
            add_msg_board(token);
        end;
    end;
    message.error=2;
    message.text='Tinker error.';
    add_msg_board('ERROR: Tinker did not run successfully.');
    % set(hwin,'Pointer','arrow');
    cd(my_dir);
    return
end;

comments=textscan(w,'%s','Delimiter','\n');
lines=comments{1};
rd_grad=false;
gradients=zeros(20000,3);
poi=0;
failed=false;
msg.error=0;
msg.text='No error.';
for k=1:length(lines),
    tline=char(lines(k));
    if ~isempty(tline),
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmpi(char(args(1)),'Cartesian') && strcmpi(char(args(2)),'Gradient') && strcmpi(char(args(3)),'Breakdown'),
            rd_grad=true;
        end;
        if rd_grad && strcmpi(char(args(1)),'Anlyt'),
            if length(args)<5,
                rd_grad=false;
                failed=true;
            else
                atnum=str2double(char(args(2)));
                if atnum>poi,
                    poi=atnum;
                end;
                gradients(atnum,1)=str2double(char(args(3)));
                gradients(atnum,2)=str2double(char(args(4)));
                gradients(atnum,3)=str2double(char(args(5)));
            end;
        end;
        if strcmpi(char(args(1)),'Total') && strcmpi(char(args(2)),'Gradient') && strcmpi(char(args(3)),'Norm'),
            rd_grad=false;
        end;
    end;
end;

if failed || poi==0,
    gradients=[];
    net_grad=[];
    msg.error=3;
    msg.text='Ill-defined gradients';
    return
end;

gradients=gradients(1:poi,:);

if isfield(correspondence,'network') && ~isempty(correspondence.network),
    net_grad=zeros(length(correspondence.network),3);
    for k=1:length(correspondence.network),
        if correspondence.network(k)>poi,
            gradients=[];
            net_grad=[];
            msg.error=3;
            msg.text='Ill-defined gradients';
            return
        end;
        net_grad(k,:)=gradients(correspondence.network(k),:);
    end;
end;