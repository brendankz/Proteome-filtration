# Proteome-filtration
Series of code that takes a species' proteome in FASTA format and filters for proteins with 3 characteristics: Size, Intrinsic Disorder, and Isoelectric Point. The purpose of this procedure is to create a filtering procedure for proteins to find similar proteins to a target based on the physical characteristics of the target to discover results that may be omitted by a BLAST sequence similarity search.

## Prediction and Filtration by Disorder and Size:

The program Metaprdict (https://github.com/idptools/metapredict) allows for predictions of disorder on a per residue basis from an entire proteome in FASTA format directly in the command line

Install with:
```
$ pip install metapredict
```

Metapredict recommends a threshold of 0.3 for significant disorder, which was used to convert the spectrum of disorder values into a binary system of predicting disorder per residue so that disorder per protein can be easily defined as a percentage.
Since disorder scores are returned on a per amino acid basis, this tool can also be used to quickly quantify the size of the protein.
The following Python code was generated to consolidate the raw data returned from Metapredict and categorize disorder by protein:

```
threshold = 0.30
with open('metapredict_data.csv', 'r') as f:
    for line in f:
        disorder = line.strip().split(',')[1:]
        amino_acids = len(disorder)
        filtered_disorder = list(filter(lambda x: float(x) > threshold, disorder))
        num_filtered = len(filtered_disorder)
        percentage = num_filtered / amino_acids * 100 if amino_acids > 0 else 0
        print(f'{line.strip().split()[0]} {amino_acids} {num_filtered} {percentage:.2f}')

```

This returns a CSV file that contains the disorder, the threshold of filtration on each amino acid, and the size of the protein.
A sample return is as follows:

GeneID Z Y X.XX%

Where column 1 shows the GeneID, column 2 shows the total amount of amino acids (Z), column 3 shows the number of amino acids that satisfied the threshold (Y), and column 4 shows what the percentage of amino acids satisfying the threshold.

From there, the CSV can be filtered based on the desired size and disorder of the target range using the awk command:

```
$ awk '$2 > A {print $1}'  metapredict_data_perprotein.csv > metapredict_size.csv
$ awk '$4 > B {print $1}' metapredict_size.csv > metapredict_size_disorder.csv
```

This returns a list of the proteins that satisfy the requirements of size and disorder given.

## **Prediction and Filtration by Iso-Electric Point**

Using the Emboss package, iso-electric point can be estimated for each protein. The Emboss package can be found here: https://anaconda.org/bioconda/emboss
Install with:
```
conda install -c bioconda emboss
```
With the Emboss package, the iso-electric point prediction tool can be used for filtration
```
$ iep proteome.fasta
```
The following Python code consolidates the raw iep return for each protein, assigning a single value and attaching the GeneID:
```
import re
with open("data.iep", "r") as f:
    data = f.read()
matches = re.findall(r'IEP of (\S+) from \d+ to \d+\nIsoelectric Point = (\d+\.\d+)', data)
for gene_id, iep in matches:
    print("{} {}".format(gene_id, iep))
```
Filter IEP by :
```
awk ‘$2>{Lower Limit} && $2<{Upper Limit} {print $1}’ iep.clean.txt > iepgenesonlyfiltered.txt
```

## Consolidating Disorder, Size, and Iso-Electric Point Filtrations
Open iepgenesonlyfiltered.txt and metapredict_size_disorder.csv and copy each list. Input these separately into https://bioinfogp.cnb.csic.es/tools/venny/ to create a Venn Diagram of each filtered list and extract the similarities to generate the final filtered list of proteins.

## **Proteome Filtration by Coiled-coil Domain**

The DeepCoil program was used for coiled-coil domain predictions (https://github.com/labstructbioinf/DeepCoil)
```
$ pip3 install deepcoil
```

Plots of each protein previously filtered for size, disorder, and iso-electric point were generated using the plot function. These plots are generated as PNG images and can be locally downloaded for analysis. The images were analyzed visually for similarity to the coiled-coil domain of the target protein in domain structure and organization.

```
deepcoil [-h] -i FILE [-out_path DIR] [-n_cpu NCPU] [--gpu] [--plot]
                [--dpi DPI]
```
