
function [model1, model2] = checkingFormulasandCharges(model1 , model2, model1MetBiGG, model1mat, model2MetBiGG, model2mat )
% 
%[model1, model2] = checkingFormulasandCharges(model1 , model2, model1MetBiGG, model1mat, model2MetBiGG, model2mat )
% 
%This function takes two metabolic models as an input as well as the
%information from BiGG(Formulas, Charges, and Atomic matrices) and outputs
%edited versions of those two models where it included 2 new fields
%"altMetCharge" and "altMetFormulas" in which it specifies alternative
%charges and formulas which help atomically and charge balance reaction
%equations 
%
% INPUT:
% model1: First model to be checked
%
% model2: second model to be checked 
%
% model1MetBiGG, model1mat : output of function getbiggChargeFormula.m for
% model 1
%
% model2MetBiGG, model2mat : output of function getbiggChargeFormula.m for
% model 2
%
% OUTPUT:
% model1, model2 : edited models 


% get unbalanced reactions as a bool column vector 
[~,~,~,imBalancedBool1,~] = checkMassChargeBalance(model1);
[~,~,~,imBalancedBool2,~] = checkMassChargeBalance(model2);

%indices of nonzero values for model1 (indices of unbalanced rxns) 
[i1,~,~] = find(imBalancedBool1);

reactionFixed1 = 0;
formulaFixed1 = 0;
chargeFixed1 = 0;
x1=0;
y1=0;
z1=0;
rxnToChange = zeros(size(model1.S,1),1);
model1.altMetCharge = model1.metCharge;
model1.altMetFormulas = cell(size(model1.S,1),1);

%going through each unbalanced reaction tocheck if it can be fixed
for q = 1:length(i1)
        rxnToChange = model1.S( :, i1(q,1));
        [x1,y1,z1, model1] = checkandChangebalance(model1, rxnToChange, model1MetBiGG, model1mat);
        reactionFixed1 = reactionFixed1+x1;
        formulaFixed1=formulaFixed1+y1;
        chargeFixed1 = chargeFixed1 +z1;
end

fprintf('\n%d Charges were edited in Model 1 ',chargeFixed1);
fprintf('\n%d Formulas were edited in Model 1 ',formulaFixed1);
fprintf('\n%d Reactions were balanced in Model 1\n',reactionFixed1);


[i2,~,~] = find(imBalancedBool2);
reactionFixed2 = 0;
formulaFixed2 = 0;
chargeFixed2 = 0;
x2=0;
y2=0;
z2=0;
rxnToChange = zeros(size(model2.S,1),1);
model2.altMetCharge = model2.metCharge;
model2.altMetFormulas = cell(size(model2.S,1),1);

for q = 1:length(i2)
        rxnToChange = model2.S( :, i2(q,1));
        [x2,y2,z2, model2] = checkandChangebalance(model2, rxnToChange, model2MetBiGG, model2mat);
        reactionFixed2 = reactionFixed2+x2;
        formulaFixed2=formulaFixed2+y2;
        chargeFixed2 = chargeFixed2 +z2;
end

fprintf('\n%d Charges were edited in Model 2 ',chargeFixed2);
fprintf('\n%d Formulas were edited in Model 2 ',formulaFixed2);
fprintf('\n%d Reactions were balanced in Model 2\n',reactionFixed2);


