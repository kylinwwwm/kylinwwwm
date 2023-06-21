%MEF_DATA = load('E:\\code\\1_test\\data\\MEF_QC_all.csv');
mESC_DATA = load('E:\\code\\1_test\\data\\mESC_data.csv');
data_mean = mean(mESC_DATA,2);
data_var = var(mESC_DATA,1,2);
data_cv2 = data_var ./ data_mean.^2;
data_fano = data_var ./data_mean;
data_sk = skewness(mESC_DATA,1,2);
data_kt = kurtosis(mESC_DATA,1,2);
data_bc = ((data_sk-1).^2 + 1)./data_kt;
gene_number_non =  struct2array(load('E:\\code\\1_test\\mESC_gene_number_non'));
gene_number = 1:size(mESC_DATA,1);
for i = 1:length(gene_number_non)
    data_mean(gene_number_non(i)) = [];
end
for i = 1:length(gene_number_non)
    data_var(gene_number_non(i)) = [];
end
for i = 1:length(gene_number_non)
    data_cv2(gene_number_non(i)) = [];
end
for i = 1:length(gene_number_non)
    data_fano(gene_number_non(i)) = [];
end
for i = 1:length(gene_number_non)
    data_sk(gene_number_non(i)) = [];
end
for i = 1:length(gene_number_non)
    data_kt(gene_number_non(i)) = [];
end
for i = 1:length(gene_number_non)
    data_bc(gene_number_non(i)) = [];
end
for j = 1:length(gene_number_non)
    gene_number(gene_number_non(j)) = [];
end
data_mean = data_mean';
data_var = data_var';
data_cv2 = data_cv2';
data_fano = data_fano';
data_sk = data_sk';
data_kt = data_kt';
data_bc = data_bc';
save mESC_feature_result.mat data_mean data_var data_cv2 data_fano data_sk data_kt data_bc gene_number