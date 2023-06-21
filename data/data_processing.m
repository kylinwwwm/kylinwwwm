c57 = csvread('E:\code\1_test\data\SS3_c57_UMIs_mESC.csv');
cast = csvread('E:\code\1_test\data\SS3_cast_UMIs_mESC.csv');
mESC_total = [c57,cast];
gene_sum = 0;
cell_sum = 0;
poor_gene_number = [];
poor_cell_number = [];
for gene_number = 1:size(mESC_total,1)
    data_mean = mean(mESC_total(gene_number,:)) ;
    for cell_number = 1:size(mESC_total,2)
        a = mESC_total(gene_number,:);
        b = a(:,cell_number);
        if b ~= 0
            flag = 1;
        else
            flag = 0;
        end
        gene_sum = gene_sum + flag;
    end
   
    if gene_sum < 40 || data_mean < 2
        poor_gene_number(end + 1) = gene_number;
    end
    gene_sum = 0;
end

for cell_number = 1:size(mESC_total,2)
    for gene_number = 1:size(mESC_total,1)
        c = mESC_total(:,cell_number);
        d = c(gene_number,:);
        if d ~= 0
            flag = 1;
        else
            flag = 0;
        end
        cell_sum = cell_sum + flag;
    end
    
    if cell_sum <= 2000
        poor_cell_number(end + 1) = cell_number;
    end
    cell_sum = 0;
end
i = poor_gene_number;
mESC_total(i,:) = [];
save('mESC_total_data','mESC_total')  
        
        
 
    