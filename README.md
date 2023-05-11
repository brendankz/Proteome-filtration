# Proteome-filtration
Series of code that takes a species' proteome in FASTA format and filters for proteins with 3 characteristics: Size, Intrinsic Disorder, and Isoelectric Point. The purpose of this procedure is to create a filtering procedure for proteins to find similar proteins to a target based on the physical characteristics of the target to discover results that would be omitted by a BLAST sequence similarity search.

# Prediction and Filtration by Disorder and Size:

The program Metaprdict (https://github.com/idptools/metapredict) allows for predictions of disorder on a per residue basis from an entire proteome in FASTA format directly in command line
Install with:
$ pip install metapredict

Metapredict recommends a threshold of 0.3 for siginifcant disorder, which was used to convert the spectrum of disorder values into a binary system of predicting disorder per residue so that disorder per protien can be easily defined as a percentage.
Since disorder scores are returned on a per amino acid basis, this tool can also be used to quickly quantify the size of the protein.
The following python code was generated to consilidate the raw data returned from metapredict and categorize disorder by protein:

threshold = 0.30
with open('metapredict_data.csv', 'r') as f:
    for line in f:
        values = line.strip().split(',')[1:]
        num_values = len(values)
        filtered_values = list(filter(lambda x: float(x) > threshold, values))
        num_filtered = len(filtered_values)
        percentage = num_filtered / num_values * 100 if num_values > 0 else 0
        print(f'Line {line.strip().split()[0]} has {num_values} values {num_filtered} of which meet the threshold ({threshold}) which is {percentage:.2f}')

This returns a csv that contains the disorder, threshold of filtration on each amino acid, and the size of the protein.
A sample return is as follows:

Line GeneID has Z values, Y of which meet the threshold (0.3), which is X.XX%

This return can be simplified in the command line to create a file that contains only GeneID, size, and disorder score with the following command:

$ cut -d ' ' -f 1,3,5,6,7,8,9,10,11,12,13,14 --complement metapredict_data_perprotein.csv > metapredict_data_clean.csv

From there, the csv can be filtered based on the desired size and disorder of the target range using the awk command:

awk '$2 > A {print $1}'  metapredict_data_perprotein.csv > metapredict_size.csv
awk '$3 > B {print $1}' metapredict_size.csv > metapredict_size_disorder.csv

This returns a list of the proteins that satisfy the requirements of size and disorder given.

# Prediction and Filtration by Iso-Electric Point

Using the Emboss package, iso-electric point can be estimated for each protein. The Emboss package can be found here: https://anaconda.org/bioconda/emboss
Install with:
conda install -c bioconda emboss


