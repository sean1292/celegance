function [model_unref,impr,t] = refine_merged_model(modelA,modelB,obj_rxn,org_code)
% [model1,impr] = refine_merged_model(model)
% refines the merged model

% INPUT:
% modelA: 1st model
% modelB: 2nd column
% obj_rxn: objective reaction, may belong to either models
% org_code: 3-4 letter organism code

% OUTPUT:
% model1: new refined model
% impr: improvements made
% t: time taken to create the improved model (model1)

% COMMENTS: Does not take genes and stoichiometry into account! Needs to be included

tic
% merge two models to create a preliminary model
model_unref = mergeTwoModels(modelA,modelB,obj_rxn);
% initialize new_model
model1 = model_unref;

% % remove any unused metabolites and remove metabolite properties that exist
impr.unused_mets = model_unref.mets(sum(any(model_unref.S,3),2)==0);
[~,iunmets] = intersect(model1.mets,impr.unused_mets);
model1.mets(iunmets) = [];
model1.S(iunmets,:) = [];
model1 = update_metabolite_properties(model1,iunmets);
fprintf('%d unused metabolties were found and removed.\n',length(iunmets));

% % find duplicate genes, uses WormBase and KEGG subsequently to look up information
% filename_local = 'E:\Downloads\c_elegans.canonical_bioproject.current\functional_descriptions.txt';
fprintf('Looking for genes in WormBase...\n');
[~,genes_not_found_wb,twice_present] = getgeneinfo_WormBase('online',model1.genes);
fprintf('Looking for genes, not found in WormBase, in KEGG...');
[~,genes_not_found] = getgeneinfo_KEGG(org_code,genes_not_found_wb);
fprintf('Finished\nCombining duplicate genes...\n');
h = waitbar(0,'Combining duplicate genes...');
steps = size(twice_present,1);
grRules = model1.grRules;
rules = model1.rules;
x = strcat({'x'},num2str([1:1:length(model1.genes)]'));
x = strrep(x,' ','');
gindex1 = zeros(size(twice_present,1),1);
gindex2 = zeros(size(twice_present,1),1);
% % first fixing refers to changing gene properties based on current rules
for i=1:size(twice_present,1)
    % % first fixing reaction-gene matrix
    gindex1(i,1) = find(strcmp(model1.genes,twice_present{i,1}));
    gindex2(i,1) = find(strcmp(model1.genes,twice_present{i,2}));
    gene_vec_ori = model1.rxnGeneMat(:,strcmp(model1.genes,twice_present{i,2}));
    gene_vec_ori = gene_vec_ori + model1.rxnGeneMat(:,strcmp(model1.genes,twice_present{i,1}));
    model1.rxnGeneMat(:,strcmp(model1.genes,twice_present{i,2})) = gene_vec_ori;
    % % first fixing grRules
    grRules = strrep(grRules,twice_present{i,1},twice_present{i,2});
    % % first fixing rules
    rules = strrep(rules,x(gindex1(i,1)),x(gindex2(i,1)));
    waitbar(i/steps);
end
% % remove duplicate gene indices from model gene list & reaction-gene matrix
model1.genes(gindex1) = []; x(gindex1) = [];
model1.rxnGeneMat(:,gindex1) = [];
close(h);
% % second fixing begins
fprintf('Finished\nUpdating rules...');
h = waitbar(0,'Updating rules...');
steps = size(twice_present,1);
rules = model1.rules;
y = strcat({'x'},num2str([1:1:length(model1.genes)]'));
y = strrep(x,' ','');
% % second fixing refers to changing only the gene rule property based on updated rules
for i=1:size(y,1)
    rules = strrep(rules,x{i,1},y{i,1});
    waitbar(i/steps);
end
close(h);
fprintf('Finished\n');
% % check if duplicate genes have been updated corectly
[~,~,del_match] = check_gene_account(model_unref,model1,twice_present);
if sum(del_match==del_match)==length(del_match)
    fprintf('Genes have updated correctly and are accounted for.\n');
else
    fprintf('Genes have not been updated correctly and are accounted for.\n');
end

% % find duplicate reactions
r = 0; % % reversibility mismatch counter
kg = 0; % % gene involvement counter
s = 0; % % stoichiometry difference counter
d = 0; % % net duplicate reaction counter
c = 1; % % query reaction counter (WHILE loop)
h = waitbar(c,'Checking reactions for duplicity and reversibility....');
while c~=length(model1.rxns)
    steps = length(model1.rxns);
    S = model1.S; % % to account for h differences..e.g. A + B -> C vs. A + B -> C + h
    S(ismember(model1.mets,{'h[c]';'h[m]';'h[n]';'h[e]'}),:) = [];
%         fprintf('Looking for a duplicate copy of %s....\n',model1.rxns{c,1});
    k = 0; % % duplicate counter
    g = 0; % % mismatch gene involvement counter
    iduprxns = [];
    iremrxns = [];
    for i=c+1:length(model1.rxns)
        if sum(any(S(:,i),2)==any(S(:,c),2))==size(S,1)
            fprintf('Another copy of %s exists at %d: %s\n',model1.rxns{c,1},i,model1.rxns{i,1});
            X1 = printRxnFormula(model1,model1.rxns{c,1},false);
            fprintf('Original: %s: %s\n',model1.rxns{c,1},X1{1,1});
            X2 = printRxnFormula(model1,model1.rxns{i,1},false);
            fprintf('Duplicate: %s: %s\n',model1.rxns{i,1},X2{1,1});
            k = k+1;
            iduprxns(k,1) = i;
            d = d+1;
            dup_rxns{d,1} = model1.rxns{c,1}; % % stores original
            dup_rxns{d,2} = model1.rxns{i,1}; % % stores duplicate
            % % check if reversibilities are similar
            if (~isempty(strfind(X1{1,1},'<=>')) && isempty(strfind(X2{1,1},'<=>'))) || (isempty(strfind(X1{1,1},'<=>')) && ~isempty(strfind(X2{1,1},'<=>')))
                fprintf('Reversibility does not match.\n');
                r = r+1;
                dup_rev{r,1} = model1.rxns{c,1}; % % stores original
                dup_rev{r,2} = model1.rxns{i,1}; % % stores duplicate
            end
            % % check if stoichiometry is same
            if sum(S(:,c)==S(:,i))~=size(S,1)
                fprintf('The stoichiometry is different for duplicate.\n');
                s = s+1;
                dup_s{s,1} = model1.rxns{c,1}; % % stores original
                dup_s{s,2} = model1.rxns{i,1}; % % stores duplicate
            end
            % % check if gene reaction involvement is same
            if sum(model1.rxnGeneMat(c,:)==model1.rxnGeneMat(i,:))==size(model1.genes,1)
                fprintf('The gene reaction association is different for duplicate.\n');
                g = g+1;
                kg = kg+1;
                iremrxns(g,1) = i;
                rem_rxns{kg,1} = model1.rxns{c,1}; % % stores original
                rem_rxns{kg,2} = model1.rxns{i,1}; % % stores duplicate
            else
                fprintf('The gene reaction association is same for duplicate.\n');
            end
        end
    end
    if g~=0
        model1.rxns(iremrxns) = [];
        model1.S(:,iremrxns) = [];
        % % remove respective reaction properties
        model1 = update_reaction_properties(model1,iremrxns);
        fprintf('Round %d removed %d reactions.\n',c,k);
    end
    c = c+1;
    waitbar(c/steps);
end
close(h);

% % identify any duplicate metabolites
[~,i] = unique(model1.mets);
if ~isempty(setdiff([1:1:length(model1.mets)]',i))
    warning('Metabolite list is not unique. There may be errors.\n');
end
% % identify if all genes have a reaction associated
if ~isempty(model1.genes(sum(model1.rxnGeneMat,1)'==0))
    warning('%d genes have no reaction associated with them.\n',length(model1.genes(sum(model1.rxnGeneMat,1)'==0)));
end

fprintf('\nSUMMARY:\n---------------------------------------------------------\n');
fprintf('%d unused metabolties were found and removed.\n',length(iunmets));
fprintf('%d duplicate genes have been found, combined and removed.\n',size(twice_present,1));
fprintf('%d genes not found in WormBase or KEGG.\n',size(genes_not_found,1));
fprintf('%d duplicate reactions were found.\n',size(dup_rxns,1));
fprintf('%d duplicate reactions could be removed.\n',size(rem_rxns,1));

impr.dup_rxns = dup_rxns;
impr.dup_rev = dup_rev;
impr.dup_s = dup_s;
impr.rem_rxns = rem_rxns;
impr.dup_genes = twice_present;
impr.notFound_genes = genes_not_found;
impr.model = model1;

t = toc;