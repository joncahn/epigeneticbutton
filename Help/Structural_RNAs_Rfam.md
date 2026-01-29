# Structural RNA Database from Rfam

This guide explains how to download structural RNA sequences from the [Rfam database](https://rfam.org/) to filter from sRNA-seq data.

For more details, see the [Rfam sequence extraction documentation](https://docs.rfam.org/en/latest/sequence-extraction.html).

## Step 1: Install Dependencies

Install easel (required for sequence extraction):

```bash
conda install easel
```

## Step 2: Download Rfam Sequences

Download all Rfam FASTA files from the database:

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/Rfam.fa.gz
```

## Step 3: Decompress

```bash
gunzip *.gz
```

## Step 4: Index the FASTA File

```bash
esl-sfetch --index Rfam.fa
```

### Troubleshooting: Duplicate Sequences

If indexing fails due to duplicate FASTA sequences, deduplicate first:

```bash
awk 'BEGIN{RS=">"; FS="\n"; ORS=""} (FNR==1){next} { name=$1; seq=$0; gsub(/(^[^\n]*|)\n/,"",seq) } !(seen[name,seq]++){ print ">" $0 }' Rfam.fa > Rfam_dedup.fa
mv Rfam_dedup.fa Rfam.fa
esl-sfetch --index Rfam.fa
```

## Step 5: Create SQL Query Files

Create `.sql` files to fetch regions of interest. Example query for rRNA:

```sql
SELECT concat(fr.rfamseq_acc,'/',seq_start,'-',seq_end)
FROM full_region fr, rfamseq rf, taxonomy tx, family f
WHERE
    rf.ncbi_id = tx.ncbi_id
    AND f.rfam_acc = fr.rfam_acc
    AND fr.rfamseq_acc = rf.rfamseq_acc
    AND tx.ncbi_id = 4577           -- NCBI ID of organism (4577 = Zea mays)
    AND f.type LIKE '%rRNA%'        -- Type of features to fetch
    AND is_significant = 1;
```

## Step 6: Create Queries for All RNA Types

Create separate `.sql` files for each structural RNA type you want to filter:
- rRNA
- tRNA
- snRNA
- snoRNA
- etc.

## Step 7: Query the Rfam Database

Execute each SQL query against the public Rfam MySQL database:

```bash
mysql -urfamro -hmysql-rfam-public.ebi.ac.uk -P4497 \
    --skip-column-names --database Rfam < query.sql > accessions.txt
```

### Combining Multiple Queries

If using multiple query files, concatenate and validate the results:

```bash
# Combine all accession files and remove duplicates
cat accessions_*.txt | sort -u > unique_accessions.txt

# Validate accessions exist in Rfam.fa
while read accession; do
    if [ $(grep -c "$accession" Rfam.fa) -gt 0 ]; then
        printf "$accession\n" >> good_accessions.txt
    fi
done < unique_accessions.txt
```

## Step 8: Extract FASTA Sequences

Extract the sequences for your validated accessions:

```bash
esl-sfetch -f /path/to/Rfam.fa /path/to/good_accessions.txt > Rfam_ncRNAs.fa
```

## Step 9: Configure the Pipeline

1. Compress the output file:
   ```bash
   gzip Rfam_ncRNAs.fa
   ```

2. Add the path to your config file:
   ```yaml
   structural_rna_fafile: "/path/to/Rfam_ncRNAs.fa.gz"
   ```
