# Andersen Lab Linkage Mapping Pipeline

## Install linkagemapping package

Install git-lfs for loading large files using the command line. ```brew install git-lfs```

```r
install_linkage_mapping <- function() {
d <- tempfile(pattern = sample(LETTERS,1))
comm <- sprintf("mkdir -p %s && cd %s && git lfs clone https://github.com/AndersenLab/linkagemapping",d, d)
system(comm)
devtools::install(paste0(d, "/linkagemapping"))
}

install_linkage_mapping()
```

## Batteries Included

This package includes all data and functions necessary to complete a mapping for the phenotype of your choice using the recombinant inbred lines from *Andersen, et al. 2015* (*G3*). Included with this package are the cross and map objects for this strain set as well a `markers.rds` file containing a lookup table for the physical positions of all markers used for mapping.

## Workflow

### Input

#### Phenotype Data

This package presents a simple workflow for massively parallel linkage mapping experiments. As input, the functions take a data frame in "long" format with the following columns (column names specified below).

##### `strain`

The strain name for each individual row.

##### `condition`

The name of the condition to which that strain was exposed during the particular experiment. For example `bleomycin` or `etoposide`. This column is optional, but plot names will be affected by an absence of a `condition` column (e.g. ".mean.TOF" vs "amsacrine.mean.TOF").

##### `trait`

The trait being measured. For example `mean.TOF` or `q75.EXT`.

##### `phenotype`

The measured quantitative value for the phenotype. For example `175.657` or `18`.

#### Cross Object

You will first need a cross object to complete the mapping. Included with this package are three cross objects to choose from:

+ `N2xCB4856cross` - N2/CB4856(Hawaiian)
+ `N2xLSJ2cross` - N2/LSJ2
+ `AF16xHK104cross` - AF16/HK104 (*C. briggsae*)
+ `N2xCB4856cross_full` -N2/CB4856 QX and ECA strains with SNPs from whole-genome sequencing

The cross objects are included in the data directory of the package. In order to access the objects in a mutable fashion, they should be read in, then assigned to another variable name as below:

```r
data(N2xCB4856cross)
cross <- N2xCB4856cross
```

### Merging the Phenotype Data Into the Cross Object

Phenotype data in the format outlined above can be merged with the cross object using the mergepheno function:

```r
# Assuming the phenotype data is in a data frame named "pheno"
phenocross <- mergepheno(cross, pheno)
```

**As of this writing (7/24/2015) the genotypes for the first two sets of RIAILs are unclear and these strains do not contribute to the mapping. Consequently, it is advised to subset to only the second set of RIAILs prior to mapping when using the N2/CB4856 cross object, as below:**

```r
# Assuming the phenotype data is in a data frame named "pheno"
phenocross <- mergepheno(cross, pheno, set = 2)
```

**As of 4/10/2017 the most recent cross object, N2xCB4856cross_full, contains SNPs from whole-genome sequence data from three sets: QX1-239, QX240-598, and ECA1-667. Select the set of strains applicable for your mapping**

### Completing the Mapping

Mapping is completed by the `fsearch` function. This function performs a mapping with forward search and returns a long-format data frame of all LOD scores across all traits in the given phenotype data.

```r
fsearch(cross)
```

In order to calculate significance cutoffs, the function needs to permute the phenotypes and remap a number of times to get either false discovery rate or genome-wide error rate. This step can take a while with 1000 iterations (the default value), so it is suggested that the number of permutations for testing large datasets locally be reduced. This will reduce the accuracy of the mappings, but will expedite exploratory analyses.

```r
fsearch(cross, permutations = 10)
```

If your computer has an NVIDEA graphics card and you have installed the gputools package, you can alternatively enable `doGPU` which will massively expedite the mapping process and allow you to complete the recommended 1000 permutations.


```r
fsearch(cross, permutations = 1000, doGPU = TRUE)
```

If you only need to map a subset of a much larger trait set, you can set a value to the `phenotype` parameter. This parameter and the subsequent trait selection follow the rules outlined by the `select` and `contains` functions from the `dplyr` package.

```r
fsearch(cross, phenotype = "bleomycin", permutations = 1000, doGPU = TRUE)
```

Set the markerset to be used for marker conversion at the end of mapping. **If you use the N2xCB4856cross_full cross object, set the markerset to NA**. Default is "N2xCB4856", other options are "N2xLSJ2" or "AF16xHK104".

```r
fsearch(cross, phenotype = "bleomycin", permutations = 1000, doGPU = TRUE, markerset = NA)
```

#### Choosing a Thresholding Method

Without overwhelming you with the details, you should generally choose false discovery rate (FDR) as your thresholding method if you are mapping more than 100 traits. Otherwise use genome-wide error rate (GWER).

Traits > 100

```r
fsearch(cross, threshold = "FDR")
```

Traits < 100

```r
fsearch(cross, phenotype = "bleomycin", threshold = "GWER")
```

### Annotating LOD scores

The final step in completing a mapping is to annotate all peak LOD scores with confidence interval bounds, effect size, and variance explained. This can be completed as below:

```r
annotate_lods(map, cross)
```

### Example

Below is a complete example of a mapping script run start to finish

```r
# Install and load the package
devtools::install_github("AndersenLab/linkagemapping")
library("linkagemapping")

# Get the cross object
data("N2xCB4856cross")
cross <- N2xCB4856cross

# Get the phenotype data
pheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/PastMembers/Tyler/ForTrip/RIAILs2_processed.rds")

# Merge the cross object and the phenotype data
cross <- mergepheno(cross, pheno)

# Perform a mapping with only 10 iterations of the phenotype data for FDR calc
map <- fsearch(cross, permutations = 10)

# Annotate the LOD scores
annotatedlods <- annotate_lods(map, cross)
```

# Marker liftover

In order too keep marker positions up to date, it will be necessary to use Dan's [liftover utilities](https://github.com/AndersenLab/liftover-utils) occasionally. For the marker positions to be updated, the markers data files in the `/data` subdirectory of the package must be updated as well as the internal markers files in `R/sysdata.rda`. Use the `use_data` function from the devtools package in order to update both the user accessible as well as the internal data sets. Additionally, make sure to update the documentation (in `R/datasets.R`) to reflect the new genome build from which the marker positions are taken.
