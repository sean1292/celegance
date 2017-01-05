function [mismatch_id,mismatch_bigg] = name_compare(name,biggid_)
% [mismatch_id,mismatch_bigg] = name_compare(name,biggid_)
% identifies metabolites which have different names or different ids;
% 1. outputs must be empty
% 2. make sure information of compartments is removed (c,n,m,e,g,y,p)

c = 1;i = 0;
mismatch_id = []; mismatch_bigg=[];
while c+1~=length(name)
    if strcmp(name{c,1},name{c+1,1})
        if strcmp(biggid_{c,1},biggid_{c+1,1})==0
            i = i+1;
            mismatch_bigg{i,1} = biggid_{c,1};
            mismatch_id(i,1) = c;
        end
    end
    c = c+1;
end