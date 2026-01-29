# Gene Ontology Database Setup

This guide explains how to create a Gene Ontology (GO) database for use with the EPICC pipeline.

## Using an Existing Database

If you already have a database from AnnotationForge, copy the file to the `genomes/<ref_genome>/GO/` directory for the analysis to work automatically.

For some species, an annotation package can be created with `AnnotationForge::makeOrgPackageFromNCBI()`. See the [AnnotationForge documentation](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html) for details.

## Creating a New Database

To create a new package de novo, two files are required:

| File | Description |
|------|-------------|
| `<ref>_infoGO.tab` | Links genes to GO terms (originally GAF format) |
| `<ref>_genes_info.tab` | Gene information (coordinates, descriptions) |

### Obtaining the GAF File

The GAF file can usually be downloaded from:
- NCBI
- Species-specific community resources (e.g., [TAIR](https://www.arabidopsis.org/) for Arabidopsis)
- [GO Consortium downloads](https://geneontology.org/docs/download-go-annotations/)

### Example: B73 (NAM/v5)

#### File 1: GO Annotations

```bash
awk '$1 !~ /^!/' maize.B73.AGPv2.aggregate.gaf.gz > data/B73_v5_infoGO.tab
```

**Required columns:**
- Column 1: Gene IDs (e.g., `AT1G00010`)
- Column 6: GO terms (e.g., `GO:00001`)
- Column 10: Evidence codes (e.g., `IEA`)

> **Note:** If your file has different column positions, edit `R_build_GO_database.R` to select the correct columns for the `fGO` table.

#### File 2: Gene Information

```bash
# Create header
printf "Chr\tStart\tEnd\tGID\tType\tDescription\n" > data/B73_v5_genes_info.tab

# Extract gene information from GFF
awk -v OFS="\t" '$3=="gene" {print $1,$4-1,$5,$9,".",$7}' genomes/B73_v5/B73_v5.gff \
  | awk -F"[:;=]" -v OFS="\t" '{print $1,$2,$4,$6}' \
  | awk -v OFS="\t" '{print $1,$2,$3,$5,$6,$7}' >> data/B73_v5_genes_info.tab
```

> **Important:** The `GID` column must match the gene IDs in File 1.

### Building the Database

Run the `R_build_GO_database.R` script with the following arguments:

```bash
script="scripts/R_build_GO_database.R"
infofile="B73_v5_infoGO.tab"      # Modified FILE1
genefile="B73_v5_genes_info.tab"  # Modified FILE2
ref_genome="B73_v5"               # Reference genome name (matches sample file)
genus="Zea"                       # Genus (capitalize first letter)
species="mays"                    # Species (lowercase)
ncbiID="4577"                     # NCBI taxonomy ID

Rscript ${script} ${infofile} ${genefile} ${ref_genome} ${genus} ${species} ${ncbiID}
```

### Configuration

Update the `GOdatabase` entry in your config file:

```yaml
GOdatabase: "org.Zmays.eg.db"
```

The naming convention is: `org.<FirstLetterGenus><species>.eg.db`

### Alternative: Using Snakemake

You can also build the database through Snakemake after filling in the config file:

```bash
snakemake --cores 1 genomes/<ref_genome>/GO/<dbname>

# Example:
snakemake --cores 1 genomes/ColCEN/GO/org.Zmays.eg.db
```

## When to Run

- **Before analysis:** Create the GO database beforehand for automatic GO analysis integration
- **After analysis:** If created later, run Snakemake with the GO file as a target to generate results
