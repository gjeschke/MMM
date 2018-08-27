function scramble = scramble_RNA_stems(stems)
% scramble = scramble_RNA_stems(stems)
%
% generates a random order for allowed preceding and following fragments of
% an RNA stem given the struct stems from an RNA library
% used for avoiding statistical bias in stem attachment
%
% see generate_RNA_stem_library
%
% G. Jeschke, 26.12.2017

scramble.previous = cell(1,length(stems));
scramble.next = cell(1,length(stems));
for ks = 1:length(stems)
    np = length(stems(ks).previous);
    scr = rand(1,np);
    [~,poivec] = sort(scr);
    scramble.previous{ks} = poivec;
    nn = length(stems(ks).next);
    scr = rand(1,nn);
    [~,poivec] = sort(scr);
    scramble.next{ks} = poivec;
end