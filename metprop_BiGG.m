function prop = metprop_BiGG(all_bigg)

options = weboptions('time',10,'RequestMethod','GET');
all_bigg = regexprep(all_bigg,'\[c\]|\[e\]|\[m\]|\[n\]','');
for i=1:length(all_bigg)
    fprintf('%d. Find in BiGG: %s...\n',i,all_bigg{i,1});
    try
        data = webread(strcat('http://bigg.ucsd.edu/api/v2/universal/metabolites/',all_bigg{i,1}),options);
    catch
        data.name = 'Not found';
        data.charges = 100000;
        data.formulae = 'Not found';
        data.database_links.KEGGCompound = 'Not found';
        data.bigg_id = 'Not found';
    end
    met_name{i,1} = data.name;
    met_charge{i,1} = data.charges;
    met_formula{i,1} = data.formulae;
    if isfield(data.database_links,'KEGGCompound')
        met_KEGG{i,1} = data.database_links.KEGGCompound;
    else
        met_KEGG{i,1} = 'Not found';
    end
    met_biggid{i,1} = data.bigg_id;
end
prop.names = met_name;
prop.charges = met_charge;
prop.formula = met_formula;
prop.KEGGID = met_KEGG;
prop.BiGG = met_biggid;
    
    