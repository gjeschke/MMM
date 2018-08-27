function anchor = get_C3p_anchor(adr)

atom_list = cell(1,12);
atom_list{1} = 'P';
atom_list{2} = 'OP1';
atom_list{3} = 'OP2';
atom_list{4} = 'O5''';
atom_list{5} = 'C5''';
atom_list{6} = 'C4''';
atom_list{7} = 'O4''';
atom_list{8} = 'C3''';
atom_list{9} = 'O3''';
atom_list{10} = 'C2''';
atom_list{11} = 'O2''';
atom_list{12} = 'C1''';

anchor = zeros(12,3);

for k = 1:length(atom_list)
    [msg,coor] = get_object(sprintf('%s.%s',adr,atom_list{k}),'coor');
    if msg.error,
        anchor = [];
        add_msg_board(sprintf('ERROR: For C3''-terminal anchor residue %s the %s atom does not exist.',adr,atom_list{k}));
        add_msg_board(msg.text);
    %     set(hfig,'Pointer','arrow');
        return
    else
        anchor(k,:) = coor;
    end;
end;