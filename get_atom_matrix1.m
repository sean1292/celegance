function atoms_mat = get_atom_matrix1(formula)
% atoms_mat = create_atom_matrix(model,mets)
% creates a cell matrix which can be evaluated using eval to determine atom balancing
% User-defined functions used:
% 1. create_atom_numbers
%
% INPUT:
% formula: formula to be translated into atom matrix
% disp_flag: true, print summary; false, don't print summary
%
% OUTPUT:
% atoms_mat: a cell matrix with atom numbers listed as string
% Columns = mets, Rows = elements
%

elements = {'As';'C';'Ca';'Cl';'Co';'Cu';'Fe';'H';'I';'K';'Mg';'Mo';'N';'Na';'O';'P';'R';'S';'Se';'Zn'};
atoms_mat = cell(1,size(elements,1));
splitFormula =regexp(formula, '([A-Z][a-z]?\d?\d?\d?)', 'tokens');
if isempty(splitFormula)~= 1
    
    for i=1:length(splitFormula)
        element=[];
        quantity=[];
        atom = regexp(splitFormula{1,i}, '([A-Z][a-z]?)(\d?\d?\d?)', 'tokens');

        if isempty(atom)~= 1
            if isempty(atom{1,1})~= 1
                if length(atom{1,1}{1,1}) == 2
                    element  = char(atom{1,1}{1,1}{1,1});
                    quantity = char(atom{1,1}{1,1}{1,2});
                else
                    element = char(atom{1,1}{1,1}{1,1});
                end
            end
        end
        
        if isempty(element)~=1
            for k=1:length(elements)
                if strcmp(element, elements{k,1})
                    if isempty(quantity)~=1
                        atoms_mat{1,k}= quantity;
                    else
                        atoms_mat{1,k}= '1';
                    end
                end
            end
        end
        
    end
    
end
