function [energy,rms_grad,grad_norm,msg] = get_Tinker_minimize(output)
% function [energy,rms_grad,grad_norm,msg] = get_Tinker_minimize(output)
%
% analyzes Tinker minimize output

msg.error=1;
msg.text='No result line';
energy=[];
rms_grad=[];
grad_norm=[];

comments=textscan(output,'%s','Delimiter','\n');
lines=comments{1};
for k=1:length(lines),
    tline=char(lines(k));
    if ~isempty(tline),
        myline = textscan(tline,'%s');
        args=myline{1};
        if strcmpi(char(args(1)),'LBFGS'),
            msg.error=0;
            msg.text=char(args(3));
        end;
        if strcmpi(char(args(1)),'Final') && strcmpi(char(args(2)),'Function') && strcmpi(char(args(3)),'Value'),
            energy=str2double(char(args(5)));
        end;
        if strcmpi(char(args(1)),'Final') && strcmpi(char(args(2)),'RMS') && strcmpi(char(args(3)),'Gradient'),
            rms_grad=str2double(char(args(5)));
        end;
        if strcmpi(char(args(1)),'Final') && strcmpi(char(args(2)),'Gradient') && strcmpi(char(args(3)),'Norm'),
            grad_norm=str2double(char(args(5)));
        end;
    end;
end;