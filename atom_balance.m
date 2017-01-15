function [rxn_bal,stat] = atom_balance(model,rxns)
% [rxn_bal,t] = atom_balance(model,rxns)
% gives the elemental imbalance in the rxns
% User-defined function used:
% 1. create_atom_numbers
% 2. create_atom_matrix
%
% INPUT:
% model: model to be checked
% rxns: list of reactions, belonging to model, to be checked
%
% OUTPUT:
% rxn_bal: elemental imbalance in reactions
% all_elements = {'C';'Ca';'Cl';'Co';'Cu';'Fe';'H';'I';'K';'Mg';'Mo';'N';'Na';'O';'P';'R';'S';'Se';'Zn'}';
% rows = rxns, columns = all_elements
% stat: balance stats (imbalance = 0, balance = 1, not checked = -1)
%
if nargin < 2
    rxns = model.rxns;
end
rxn_bal = zeros(size(model.rxns,1),1);
fprintf('Rxn no. \tReaction name\t\t\tComment\n------------------------------------------------\n');
sym('n');
for i=1:length(rxns)
    rxn_index = find(strcmp(model.rxns,rxns{i,1}));
    met_index = find(model.S(:,rxn_index)~=0);
    mets = model.mets(met_index);
    atoms_mat = create_atom_matrix(model,mets);
%     printRxnFormula(model,model.rxns(i));
    if sum(cellfun(@isempty,atoms_mat),2)==size(atoms_mat,2)
        rxn_bal(i,1) = -1;
        fprintf('%d. %s: One or more metabolites missing formula.\n',i,rxns{i,1});
    else
        Sj = zeros(1,size(atoms_mat,2));
        n=50;
        for j=1:length(met_index)
            atoms_mat(j,cellfun(@isempty,atoms_mat(j,:))) = {'0'};
            atoms_mat(j,strcmp(atoms_mat(j,:),'*n')) = {'0'};
            Sj = Sj + model.S(met_index(j),rxn_index)*cellfun(@eval,atoms_mat(j,:));
        end
        stat(i,:) = Sj;
        if sum(stat(i,:)) == 0
            rxn_bal(i,1) = 1;
            fprintf('%d. %s: Atoms balanced.\n',i,rxns{i,1});
        else
            rxn_bal(i,1) = 0;
            fprintf('%d. %s: Atoms not balanced.\n',i,rxns{i,1});
        end
    end
end