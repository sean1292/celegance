function [commonGeneRxnsModel1, commonGeneRxnsModel2, rxnsModel1, rxnsModel2] = rxnsFromGenes(model1, model2)
    genesModel1Not2 = setdiff(model1.genes, model2.genes);
    genesModel2not1 = setdiff(model2.genes, model1.genes);
    genesPresentInBoth = intersect(model1.genes, model2.genes);

    %reactions assosiated with same genes Kaleta
    rxnsBoth2S = findRxnsFromGenes(model2, genesPresentInBoth);

    %reactions assosiated with same genes Kaleta
    rxnsBoth1S = findRxnsFromGenes(model1, genesPresentInBoth);

    %reactions unique to each model according to their genes
    rxns1 = findRxnsFromGenes(model1, genes1Not2);


    rxns2 = findRxnsFromGenes(model2, genes2Not1);

    rxnsBoth2Finished = fieldnames(rxnsBoth2S)
    for  i = 1:numel(rxnsBoth2Finished)
        gene = rxnsBoth2S.(rxnsBoth2Finished{i});
        str='';
        for j = 1:size(gene,1)
            rxns(1,j) = {gene{j,1}};
        end 
        rxns = [rxns{:}];
        str = strjoin(rxns, ', ');
        rxnsBoth2Finished(i,2) = {str};
        rxns(1,:)=[];
    end

    rxnsBoth1Finished = fieldnames(rxnsBoth1S)
    for  i = 1:numel(rxnsBoth1Finished)
        gene = rxnsBoth1S.(rxnsBoth1Finished{i});
        str='';
        for j = 1:size(gene,1)
            rxns(1,j) = {gene{j,1}};
        end 
        rxns = [rxns{:}];
        str = strjoin(rxns, ', ');
        rxnsBoth1Finished(i,2) = {str};
        rxns(1,:)=[];
    end

    rxns1Finished = fieldnames(rxns1)
    for  i = 1:numel(rxns1Finished)
        gene = rxns1.(rxns1Finished{i});
        str='';
        for j = 1:size(gene,1)
            rxns(1,j) = {gene{j,1}};
        end 
        rxns = [rxns{:}];
        str = strjoin(rxns, ', ');
        rxns1Finished(i,2) = {str};
        rxns(1,:)=[];
    end

    rxns2Finished = fieldnames(rxns2)
    for  i = 1:numel(rxns2Finished)
        gene = rxns2.(rxns2Finished{i});
        str='';
        for j = 1:size(gene,1)
            rxns(1,j) = {gene{j,1}};
        end 
        rxns = [rxns{:}];
        str = strjoin(rxns, ', ');
        rxns2Finished(i,2) = {str};
        rxns(1,:)=[];
end



