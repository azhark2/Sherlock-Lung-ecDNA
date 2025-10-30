drivers = pd.read_csv("~/iCloud/dev/Sherlock-Lung/results2/sherlock_driver_mutations.tsv", sep="\t")
scna = pd.read_csv("~/iCloud/dev/Sherlock-Lung/results2/sherlock_focal_scna.tsv", sep="\t")
drivers["Gene"] = drivers["Hugo_Symbol"]


deletions = scna[scna['Alteration'] == 'Del']
amps = scna[scna['Alteration'] == 'Amp']

gene_set = set(deletions["Gene"]) | set(drivers["Gene"])

result_dict = {}
for gene in gene_set: #for every gene that was mutated/deleted

    #first iterate through SCNA deletions
    scna2 = deletions[deletions["Gene"] == gene]
    for i, row in scna2.iterrows():
        key = (row["Tumor_Barcode"], gene)
        if key not in result_dict:
            result_dict[key] = 0
        result_dict[key] += 1
    
    #then iterate through driver mutations
    drivers2 = drivers[drivers["Gene"] == gene]
    for i, row in drivers2.iterrows():
        key = (row["Tumor_Barcode"], gene)
        #check that the gene wasn't amplified 
        
        if key not in result_dict:
            result_dict[key] = 0
        result_dict[key] += 1

s = []
g = []
f = []
for k,v in result_dict.items():
    # if v == 2:
    #     print(k,v)
    s.append(k[0])
    g.append(k[1])
    f.append(v)

df = pd.DataFrame({"Tumor_Barcode": s, "Gene": g, "Frequency": f})
df["Frequency"].value_counts()

#now integrate ecDNA status
ecDNA = pd.read_csv("~/iCloud/dev/Sherlock-Lung/results2/ecDNA-annotations.txt", sep="\t")
df = pd.merge(df, ecDNA, on="Tumor_Barcode")

for g in df["Gene"].unique():
    df2 = df[df["Gene"] == g]
    num_duplicates = df2["Tumor_Barcode"].duplicated().sum()
    assert(num_duplicates == 0) 


df.to_csv("~/iCloud/dev/Sherlock-Lung/results2/drivers-scna-merged2.tsv", sep="\t", index=False)