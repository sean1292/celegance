function [x, formulasChanged, chargesChanged, model] = checkandChangebalance(model, reactionMets, biggList, biggMat)
% INPUT:
% model: metabolic model for which we would like to check balancing
% reactionMets: reaction column vector that is going to checked 
% OUTPUT:
% x:  1 or 0 depending on whether the reaction was changed (balanced) by
% the algorithm

x = 0;
formulasChanged = 0; 
chargesChanged = 0; 
[reactants, ~, number] = find(reactionMets < 0);
[products, ~, number1] = find(reactionMets > 0);

if isempty(reactants) || isempty(products)
    return
end
reactantBigg = cell(size(reactants,1),size(biggList,2));
reactantMats = zeros(size(reactants,1),size(biggMat,2));
for i = 1:size(reactants,1)
reactantBigg(i,:) = biggList(reactants(i,1), :);
end
for i = 1:size(reactants,1)
reactantMats(i,:)= biggMat(reactants(i,1), :); 
end

%[reactantBigg, reactantMats] = getbiggChargeFormula(model.mets(reactants),0);

%checking whether all of the metabolites had charge and 
for i =1:size(reactantBigg, 1)
   if strcmp('Not found', reactantBigg{i,3}) || strcmp('Not found', reactantBigg{i,2}) || strcmp('Multiple Formulas Not Resolevd', reactantBigg{i,3})
        return;
   end  
end
for i = 1:size(reactantBigg, 1) 
	(reactantMats(i,:)).*(number(i,1));
end
totReactants = sum(reactantMats);

productBigg = cell(size(products,1),size(biggList,2));
productMats = zeros(size(products,1),size(biggMat,2));

for i = 1:size(products,1)
productBigg(i,:) = biggList(products(i,1), :);
end

for i = 1:size(products,1)
productMats(i,:)= biggMat(products(i,1), :); 
end

for i =1:size(productBigg, 1)
   if strcmp('Not found', productBigg{i,3}) || strcmp('Not found', productBigg{i,2})|| strcmp('Multiple Formulas Not Resolevd', productBigg{i,3})
        return;
   end  
end

for i = 1:size(productBigg, 1) 
	productMats(i,:).*(number1(i,1));
end
totCh1 =0;
totCh2 =0;

for i=1:size(reactantBigg,1)
    totCh1 = totCh1 + str2double(reactantBigg{i,2})*(number(i,1));
end

for i=1:size(productBigg,1)
    totCh2 = totCh2 + str2double(productBigg{i,2})*(number1(i,1));
end

totProducts = sum(productMats);


if isequal(totProducts, totReactants) && totCh1==totCh2
    x = 1;
    for k=1:size(reactants, 1)
       if strcmp(string(model.metFormulas(reactants(k,1), 1)),reactantBigg{k,3})~=1
            formulasChanged = formulasChanged + 1;
            model.altMetFormulas{reactants(k,1), 1} = reactantBigg{k,3};
       end
       if model.metCharge(reactants(k,1), 1) ~= int32(str2double(reactantBigg{k,2}))
            chargesChanged = chargesChanged + 1;
            model.altMetCharge(reactants(k,1), 1) = int32(str2double(reactantBigg{k,2}));
       end
       
    end
    for k=1:size(products, 1)
       if strcmp(string(model.metFormulas(products(k,1), 1)),productBigg{k,3})~=1
            formulasChanged = formulasChanged + 1; 
            model.altMetFormulas{products(k,1), 1} = productBigg{k,3};
       end
       if model.metCharge(products(k,1), 1) ~= int32(str2double(productBigg{k,2}))
            chargesChanged = chargesChanged + 1;
            model.altMetCharge(products(k,1), 1) = int32(str2double(productBigg{k,2}));
       end
    end
end 


