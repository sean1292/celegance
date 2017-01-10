function gene_list = getgeneinfo_KEGG(org_code,gene_entry,web_options)

if nargin < 3
    web_options = weboptions('time',10);
end
gene_list = cell(length(gene_entry),1);
for i=1:length(gene_entry)
    A = webread(strcat('http://rest.kegg.jp/find/',org_code,'/',gene_entry{i,1}),web_options);
    A = getallgenes_KEGG(A,org_code,false);
    A = A(strcmp(A(:,1),gene_entry{i,1}),:);
    gene_list{i,1} = A;
    fprintf('%d.Query: %s....\n',i,gene_entry{i,1});
end