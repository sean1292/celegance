function [metBiGG, atom_mat] = getbiggChargeFormula(all_biggid, disp)
% metBiGG = getbiggChargeFormula(all_bigg, 1)
% INPUT:
% all_bigg: list of Bigg compound ids
% disp: flag - 1 to show outpput - 0 - not to show
% OUTPUT:
% metBiGG: list of Chemical Formulas and Charges from
%          Bigg
% atom_mat - matrix containing #s of each chemical element
%            in each fomula retrieved from Bigg
tic;
metBiGG = cell(length(all_biggid),9);
options = weboptions('Timeout',30, 'RequestMethod' , 'GET');

notInbigg = 0;
noCharge=0;
noFormula=0;
inchiFound = 0;
chResolved = 0;
nochResolved = 0;
formResolved = 0;
noformResolved = 0;
gibbsFound= 0;
found = 0;
atomMat= 0;
noAtomMat= 0;
notSameFormulaMetacyc = 0;

%check all Bigg IDs
for i=1:length(all_biggid)
    metBiGG{i,1}= all_biggid{i,1};
    inchi =[];
    
    try
        %if not same id as previous one
        
        data= webread(strcat('http://bigg.ucsd.edu/api/v2/universal/metabolites/', all_biggid{i,1}(1:end-3)),options);
        
        
        
        %if there is a charge in BIGG
        if isempty(data.charges)~=1
            
            %getting all charges
            for k=1:length(data.charges)
                if k == length(data.charges) && length(data.charges)==1
                    metBiGG{i,2} = strcat(metBiGG{i,2}, num2str(data.charges(k,1)));
                else
                    metBiGG{i,2} = strcat(metBiGG{i,2}, num2str(data.charges(k,1)), ', ');
                    metBiGG{i,4} = length(data.charges); %#of charges present
                    if k == length(data.charges)
                        if isfield(data.database_links,'MetaNetX_MNX_Chemical')
                            metaNetX = webread( data.database_links.MetaNetX_MNX_Chemical.link);
                            regex = regexp(metaNetX, 'InChI=(.+)</td></tr>   <tr><td', 'tokens');
                            if isempty(regex)~=1
                                inchi = regex{1};
                            end
                        end
                        if isempty(inchi)~=1
                            metBiGG{i,6} = inchi;
                            inchiFound = inchiFound +1;
                        end

                        if  isempty(metBiGG{i,6})~=1
                            metBiGG{i,7} = metBiGG{i,2};
                            metBiGG{i,2} = getChargeFromInChI(char(metBiGG{i,6}));
                            chResolved = chResolved+1; %when multiple charges present
                        end
                    end
                end
            end
        else
            noCharge = noCharge+1;
            if isfield(data.database_links,'MetaNetX_MNX_Chemical')
                metaNetX = webread( data.database_links.MetaNetX_MNX_Chemical.link);
                regex = regexp(metaNetX, 'InChI=(.+)</td></tr>   <tr><td', 'tokens');
                if isempty(regex)~=1
                    inchi = regex{1};
                end
            end
            if isempty(inchi)~=1
                metBiGG{i,6} = inchi;
                inchiFound = inchiFound +1;
            end
            
            
            metBiGG{i,2} = 'No Charge In BiGG';
            if  isempty(metBiGG{i,6})~=1
                metBiGG{i,7} = metBiGG{i,2};
                metBiGG{i,2} = getChargeFromInChI(char(metBiGG{i,6}));
                nochResolved = nochResolved+1;
            end
        end
        
        %if there is a Formula in BIGG
        if isempty(data.formulae)~=1
            
            %getting all charges
            for x=1:length(data.formulae)
                if x == length(data.formulae) && length(data.formulae)==1
                    metBiGG{i,3} = strcat(metBiGG{i,3},data.formulae{x,1});
                else
                    metBiGG{i,3} = strcat(metBiGG{i,3}, data.formulae{x,1}, ', ');
                    metBiGG{i,5} = length(data.formulae);
                    if x == length(data.formulae)
                        if isfield(data.database_links,'MetaNetX_MNX_Chemical')
                            metaNetX = webread( data.database_links.MetaNetX_MNX_Chemical.link);
                            regex = regexp(metaNetX, 'InChI=(.+)</td></tr>   <tr><td', 'tokens');
                            if isempty(regex)~=1
                                inchi = regex{1};
                            end
                        end
                        if isempty(inchi)~=1
                            metBiGG{i,6} = inchi;
                            inchiFound = inchiFound +1;
                        end
                        %when more than one formula present
                        if  isempty(metBiGG{i,6})~=1
                            metBiGG{i,8} = metBiGG{i,3};
                            metBiGG{i,3}=getFormulaFromInChI(char(metBiGG{i,6}));
                            formResolved = formResolved+1;
                        else
                            metBiGG{i,8} = metBiGG{i,3};
                            metBiGG{i,3} = 'Multiple Formulas Not Resolevd';
                        end
                    end
                end
            end
        %trying to get formula for metabolites that are not in bigg    
        else
            noFormula = noFormula+1;
            if isfield(data.database_links,'MetaNetX_MNX_Chemical')
                metaNetX = webread( data.database_links.MetaNetX_MNX_Chemical.link);
                regex = regexp(metaNetX, 'InChI=(.+)</td></tr>   <tr><td', 'tokens');
                if isempty(regex)~=1
                    inchi = regex{1};
                end
            end
            if isempty(inchi)~=1
                metBiGG{i,6} = inchi;
                inchiFound = inchiFound +1;
            end
            
            
            metBiGG{i,8} = 'No formula In BiGG';
            
            if  isempty(metBiGG{i,6})~=1
                metBiGG{i,3} = getFormulaFromInChI(char(metBiGG{i,6}));
                noformResolved = noformResolved+1;
            end
        end
        
        %checking for Gibbs Free energy while checking if same chemical formula 
        if isfield(data.database_links,'BioCyc')
            gibbs{1}= 0;
            for j=1:length(data.database_links.BioCyc)
                metaCyc = webread( strcat('https://websvc.biocyc.org/getxml?id=', data.database_links.BioCyc(j,1).id), options);
                regexGibbs = regexp(metaCyc, '<gibbs-0 datatype=.+units.+>(-?\d+\.\d+)\n.*<evidence>', 'tokens');
                regexFormula = regexp(metaCyc, '<formula concise=.([\w\s\d]+).\/>', 'tokens');
                if isempty(regexFormula)~=1
                    regexFormula = char(regexFormula{1,1});
                    regexFormula = regexFormula(regexFormula~=' ');
                end
                metacycFormulaMat = get_atom_matrix1(regexFormula);
                biggFomulaMat = get_atom_matrix1(metBiGG{i,3});
                if isempty(regexGibbs)~=1 && isequal(metacycFormulaMat, biggFomulaMat)
                    gibbs = regexGibbs{1};
                    gibbsFound = gibbsFound+1;
                    metBiGG{i,9}= gibbs;
                elseif isequal(metacycFormulaMat, biggFomulaMat)~=1

                    notSameFormulaMetacyc =notSameFormulaMetacyc+1;
                end
            end
        end
        
        %checking if display flag is set to ON
        if disp == 1
            fprintf('\n%d. Looking for Formula and Charge for %s: ',i,all_biggid{i,1}(1:end-3));
            if isempty(data.formulae)~=1
                fprintf('\nFormula: ');
                for y=1:length(data.formulae)
                    fprintf('%s ',data.formulae{y,1});
                    if y==length(data.formulae)
                        fprintf('\n');
                    end
                end
            end
            if isempty(data.charges)~=1
                fprintf('\nCharge(s): ');
                for z=1:length(data.charges)
                    fprintf('%s ',num2str(data.charges(z,1)));
                    if z==length(data.charges)
                        fprintf('\n');
                    end
                end
            end
            if isempty(gibbs)~=1
                fprintf('Gibbs Fee Energy: %s', gibbs{1});
            end
        end
        
        
        
        %case where webread returns error
    catch
        if disp == 1
            notInbigg = notInbigg+1;
            metBiGG{i,2}= 'Not found';
            metBiGG{i,3}= 'Not found';
            fprintf('\n%d. Looking for Formula and Charge for %s: ',i,all_biggid{i,1}(1:end-3));
            fprintf('\n%s.', 'Not found')
        end
    end
    if i == length(all_biggid)
        fprintf('\nSummary:\n');
        fprintf('%d Metabolites were not in the BiGG database\n',notInbigg);
        fprintf('%d Metabolites did not have charge specified in the BiGG database (%d Resolved using InChi)\n',noCharge, nochResolved);
        fprintf('%d Metabolites did not have a formula specified in the BiGG database (%d Resolved using InChi)\n',noFormula, noformResolved);
        fprintf('%d Charges Resolved using InChi ID (where more than one charge present in BiGG)\n',chResolved);
        fprintf('%d Formulas Resolved using InChi ID (where more than one formula present in BiGG)\n',formResolved);
        fprintf('Found %d Gibbs Free Energy Values(%d did not have the same formula in Metacyc)\n',gibbsFound, notSameFormulaMetacyc);
    end
end
empty= cell(1,20);

for i=1:size(metBiGG,1)
    if isempty(metBiGG{i,3})~=1
        if strcmp('Not found', metBiGG{i,3})~=1 && (strcmp('Multiple Formulas Not Resolevd', metBiGG{i,3})~=1)
            if i == 1
                atomMat= atomMat+1;
                atom_mat = get_atom_matrix1(metBiGG{i,3});
            else
                atomMat= atomMat+1;
                atom_mat = [atom_mat; get_atom_matrix1(metBiGG{i,3})];
            end
        else
            noAtomMat= noAtomMat+1;
            if i == 1
                atom_mat = empty;
            else
                atom_mat = [atom_mat; empty];
            end
            
        end
    else
        noAtomMat= noAtomMat+1;
        if i == 1
            atom_mat = empty;
        else
            atom_mat = [atom_mat; empty];
        end
        
    end
    if i == length(metBiGG)
        fprintf(' %d Metabolites read and converted to matrix form (%d did not have a formula)\n', atomMat, noAtomMat);
    end
end
toc;