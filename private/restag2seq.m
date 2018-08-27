function seq = restag2seq(restags)

global residue_defs

seq = '';
found = true;
res = 0;
while found
    res = res + 1;
    tag = id2tag(res,restags);
    if isempty(tag)
        found = false;
    else
        aa = tag2id(upper(tag),upper(residue_defs.restags));
        if ~isempty(aa)
            seq = [seq residue_defs.single_letter_code(aa)];
        else
            seq = [seq '?'];
        end;
    end
end;