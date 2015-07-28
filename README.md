# Andersen Lab Linkage Mapping Pipeline

## Batteries Included

This package includes all data and functions necessary to complete a mapping for the phenotype of your choice using the recombinant inbred lines from *Andersen, et al. 2015* (*G3*). Included with this package are the cross and map objects for this strain set as well a `markers.rds` file containing a lookup table for the physical positions of all markers used for mapping.

## Workflow

### Input

#### Phenotype Data

This package presents a simple workflow for massively parallel linkage mapping experiments. As input, the functions take a data frame in "long" format with the following columns (column names specified below).

##### `strain`

The strain name for each individual row.

##### `condition`

The name of the condition to which that strain was exposed during the particular experiment. For example `bleomycin` or `etoposide`.

##### `trait`

The trait being measured. For example `mean.TOF` or `q75.EXT`.

##### `phenotype`

The measured quantitative value for the phenotype. For example `175.657` or `18`.

#### Cross Object

The cross object is included the data directory of the package. In order to access the object in a mutable fashion, it should be read in, then assigned to another variable name as below:

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

As of this writing (7/24/2015) the genotypes for the first two sets of RIAILs are unclear and these strains do not contribute to the mapping. Consequently, it is advised to subset to only the third set of RIAILs prior to mapping, as below:

```r
# Assuming the phenotype data is in a data frame named "pheno"
phenocross <- mergepheno(cross, pheno, set = 3)
```

### Completing the Mapping



