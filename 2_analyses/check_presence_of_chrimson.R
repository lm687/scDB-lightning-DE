all_samples_data <- list(FEB2019_G1_rep1.data,
                    FEB2019_G1_rep2.data,
                    FEB2019_G2_rep1.data,
                    FEB2019_G2_rep2.data,
                    FEB2019_G3_rep1.data,
                    FEB2019_G3_rep2.data,
                    FEB2019_G4_rep1.data,
                    FEB2019_G4_rep2.data)


gene_names_replace[chrimsonidx]


table(is.na(gene_names_replace))
# is.na(gene_names_replace[17897])
length(gene_names_replace)


rownames(all_samples$G1_rep1)[chrimsonidx]

all_samples[[1]]@assays$RNA['Chrimson']

View(data.frame(gene_names_replace))

chrimsonidx <- (which(rownames(FEB2019_G1_rep1.data) == 'FBto0000555'))
chrimsonidx
name_conversion_file[grepl('Chrimson', name_conversion_file$V2),]

find_chrimson <- function(i, name='Chrimson'){
  i[grepl(name, i)]
}

View(data.frame(rownames(FEB2019_G1_rep1@assays$RNA)))

find_chrimson(rownames(FEB2019_G1_rep1@assays$RNA))
find_chrimson(name_conversion_file$V2) ## chrimson is here
find_chrimson(rownames(FEB2019_G1_rep1.data))
FEB2019_G1_rep1@assays$RNA@counts@Dimnames[grepl('FBto0000555', FEB2019_G1_rep1@assays$RNA@counts@Dimnames)]

View(data.frame(rownames(FEB2019_G1_rep1.data)))

find_chrimson(rownames(FEB2019_G1_rep1.data), 'FBto0000555')
find_chrimson((FEB2019_G1_rep1.data@Dimnames[[1]]), 'FBto0000555')
find_chrimson((FEB2019_G1_rep1.data@Dimnames[[1]]), 'FBto0000555')
find_chrimson(dimnames(FEB2019_G1_rep1)[[1]])

table(FEB2019_G1_rep1.data['FBto0000555',]) ### always zero for chrimson counts
table(FEB2019_G1_rep2.data['FBto0000555',]) ### always zero for chrimson counts

sapply(all_samples_data, function(i) mean(i['FBto0000555',]) ) ## there isn't any chrimson anywhere
sapply(all_samples_data, function(i) mean(i['FBgn0040372',]) ) ## some random gene


gene_names_replace <- name_conversion_file[,2][match(FEB2019_G1_rep1@assays$RNA@counts@Dimnames[[1]], name_conversion_file[,1])]
gene_names_replace[grepl('Chrimson', gene_names_replace)]
