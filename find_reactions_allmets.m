function [rxns] = find_reactions_allmets(model,all_mets)

rxns = '';
for i=1:length(all_mets)
    imets = ismember(model.mets,all_mets{i,1});
    rxn_list = model.rxns(any(model.S(imets,:),1));
    if length(all_mets)==1 || i==1
        rxns = rxn_list;
    else
        rxns = intersect(rxn_list,rxns);
    end
end
