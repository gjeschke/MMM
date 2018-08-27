global general 

builddocsearchdb(general.help_files);

pfname=which('preferences.mat');

load preferences

virgin=true;

save(pfname,'virgin','user_preferences');

load third_party_references.mat

m=length(references);
auto_ref_tags='|';
for k=1:m,
    auto_ref_tags=[auto_ref_tags references(k).short '|'];
end;

save third_party_references references auto_ref_tags

