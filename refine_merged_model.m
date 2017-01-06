function [model1,impr] = refine_merged_model(model)
% [model1,impr] = refine_merged_model(model)
% refines the merged model

% INPUT: 
% model: combined model obtained from following command
% model = mergeTwoModels(eleg,icel,'BIO0100');

% OUTPUT:
% model1: new refined model
% impr: improvements made

% COMMENTS: Does not take genes into account! Needs to be included

% initialize new_model
model1 = model;

% remove any unused metabolites and remove metabolite properties that exist
impr.unused_mets = model.mets(sum(any(model.S,3),2)==0);
[~,~,iunmets] = intersect(model1.mets,impr.unused_mets);
model1.mets(iunmets) = [];
model1.S(iunmets,:) = [];
% metabolite properties
if isfield(model1,'metNames') % metabolite name
    model1.metNames(iunmets) = [];
end
if isfield(model1,'metFormulas') % metabolite formula
    model1.metFormulas(iunmets) = [];
end
if isfield(model1,'metCharge') % metabolite charge
    model1.metCharge(iunmets) = [];
end
if isfield(model1,'b') % metabolite accumulation
    model1.b(iunmets) = [];
end
if isfield(model1,'metKEGGID') % metabolite KEGG
    model1.metKEGGID(iunmets) = [];
end
if isfield(model1,'metInChIString') % metabolite InChIString
    model1.metInChIString(iunmets) = [];
end
if isfield(model1,'metChEBIID') % metabolite ChEBI
    model1.metChEBIID(iunmets) = [];
end
if isfield(model1,'metPubChemID') % metabolite PubChem
    model1.metPubChemID(iunmets) = [];
end
fprintf('%d unused metabolties were found and removed.\n',length(iunmets));

% find duplicate reactions
c = 1;
while c~=length(model1.rxns)
    conmets_c = model1.mets(model1.S(:,c)<0);
    promets_c = model1.mets(model1.S(:,c)>0);
    fprintf('Looking for a duplicate copy of %s....\n',model1.rxns{c,1});
    k = 0; % duplicate counter
    iduprxns = [];
%     fprintf('Next loop length: %d.\n',length([c+1:1:length(model1.rxns)]));
    for i=c+1:length(model1.rxns)
%         fprintf('Query: %d\tList: %d\n',c,i);
        conmets_i = model1.mets(model1.S(:,i)<0);
        promets_i = model1.mets(model1.S(:,i)>0);
        cons_int = intersect(conmets_i,conmets_c); cons_un = union(conmets_i,conmets_c);
        prod_int = intersect(promets_i,promets_c); prod_un = union(promets_i,promets_c);
        if length(cons_un)==length(cons_int) && length(prod_un)==length(prod_int) && ~isempty(cons_un)
            fprintf('Another copy of %s exists at %d: %s\n',model1.rxns{c,1},i,model1.rxns{i,1});
            X = printRxnFormula(model,model1.rxns{c,1},false);
            fprintf('Original: %s: %s\n',model1.rxns{c,1},X{1,1});
            X = printRxnFormula(model,model1.rxns{i,1},false);
            fprintf('Duplicate: %s: %s\n',model1.rxns{i,1},X{1,1});
            k = k+1;
            iduprxns(k,1) = i;
        end
    end
    if k~=0
        model1.rxns(iduprxns) = [];
        model1.S(:,iduprxns) = [];
        % reaction properties
        if isfield(model1,'rxnNames') % reaction names
            model1.rxnNames(iduprxns) = [];
        end
        if isfield(model1,'lb')
            model1.lb(iduprxns) = []; % reaction lower bounds
            model1.ub(iduprxns) = []; % reaction upper bounds
        end
        if isfield(model1,'subSystems') % reaction subsystem
            model1.subSystems(iduprxns) = [];
        end
        if isfield(model1,'rxnReferences') % reaction references
            model1.rxnReferences(iduprxns) = [];
        end
        if isfield(model1,'rxnECNumbers') % reaction EC numbers
            model1.rxnECNumbers(iduprxns) = [];
        end
        if isfield(model1,'rxnNotes') % reaction notes
            model1.rxnNotes(iduprxns) = [];
        end
        if isfield(model1,'c') % reaction objective coefficient
            model1.c(iduprxns) = [];
        end
        if isfield(model1,'rev') % reaction reversibility
            model1.rev(iduprxns) = [];
        end
        if isfield(model1,'rxnGeneMat') % reaction gene association
            model1.rxnGeneMat(iduprxns) = [];
        end
        if isfield(model1,'grRules') % reaction gene rules
            model1.grRules(iduprxns) = [];
        end
        if isfield(model1,'rules') % reaction rules
            model1.rules(iduprxns) = [];
        end
    end
    fprintf('Round %d removed %d reactions.\n',c,k);
    c = c+1;
end
impr.dup_rxns = setdiff(model.rxns,model1.rxns);
fprintf('%d duplicate reactions were removed.\n',length(impr.dup_rxns));