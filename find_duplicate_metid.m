function [dup_id,dup_name] = find_duplicate_metid(metid,name)
% [dup_id,dup_name] = find_duplicate_metid(metid,name)
% identifies duplicate metabolite ids

[~,i] = unique(metid);
i = setdiff([1:1:length(metid)],i);
dup_id = metid(i);
dup_name = name(i);

