address='[NhaA](A)<E.A>"Leu"';
%address='[1](A)"Gly"';
% address='[1](A)<L.2-3>';
%address='[1](A)<E.B>2-4.CA';
%address='[T4L]"Asn".:!';
% address='[1]|translocation<H.IVa>|';
% address='[NhaA](A)<L.2-3>';
% address='[2W8A](chain_1)"MSE"';
address='[NhaA]|bundle|';
address='[NhaA]|bundle<H.III>|';
address='[NhaA]|upper<H.IVa>|';
address='[2BHW]|Nterminal|';
address='[NhaA](A)38';

tic,
[indices,message]=resolve_address(address);
toc,


disp(message);
disp(indices);

