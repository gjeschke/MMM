clear all

initialize_MMM

global model

tic,
[message,snum]=add_pdb('1ZCD');
s_tag=sprintf('[%i]',snum);
synonym(s_tag,'NhaA');
synonym(s_tag,'EPR');
[message,snum]=add_pdb('2W8A');
synonym('[2W8A](A)','chain_1');
synonym('[2W8A](B)','chain_2');
synonym('[2W8A](C)','chain_3');
[message,snum]=add_pdb('2CUU');
synonym('[2CUU]','T4L');
synonym('[2CUU]','Hubbell');
[message,snum]=add_pdb('CuGMP_30');
s_tag=sprintf('[%i]',snum);
synonym(s_tag,'CuGMP');
message=add_pdb('2K98');
synonym('[2K98]','NMR');
message=add_pdb('1QJP');
message=add_pdb('2BHW');
synonym('[2BHW]','2BHW');
synonym('[2BHW]','LHCIIb');
synonym('[2BHW]','LHCII');
toc,

message=secondary('[NhaA](:)','helix','I',[12,30]);
message=secondary('[NhaA](:)','helix','II',[35,41]);
message=secondary('[NhaA](:)','helix','III',[59,85]);
message=secondary('[NhaA](:)','helix','IVa',[95,103]);
message=secondary('[NhaA](:)','helix','IVb',[107,116]);
message=secondary('[NhaA](:)','helix','IV',[95,116]);
message=secondary('[NhaA](:)','sheet','A',[44,50]);
message=secondary('[NhaA](:)','sheet','B',[53,58]);
message=secondary('[NhaA](:)','loop','II-III',[31,34]);

message=domain('bundle',1,'[NhaA](A)<H.I>');
message=domain('bundle',0,'[NhaA](A)<H.II>');
message=domain('bundle',0,'[NhaA](A)<H.III>');
message=domain('bundle',0,'[NhaA](A)<H.IV>');
message=domain('upper',1,'[NhaA](A)<H.IVa>');
message=domain('Nterminal',1,'[2BHW]9-58');

