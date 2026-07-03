# Example Data - Drug-Resistant Bacteria Antigen Proteins

This directory contains FASTA files with known antigen and virulence protein sequences from clinically important drug-resistant bacteria. These sequences are intended for testing and demonstration of the BacterialAntigenFinder pipeline.

## Data Files

### 1. `pseudomonas_aeruginosa_antigens.fasta`

Antigen and virulence proteins from **Pseudomonas aeruginosa** (铜绿假单胞菌), a Gram-negative opportunistic pathogen notorious for multidrug resistance, especially in hospital-acquired infections.

| Protein ID | Description | Length (aa) | UniProt Accession | Function |
|---|---|---|---|---|
| OprF_PSEAE | Outer membrane protein F | 350 | P13794 | Major porin; maintains membrane integrity and virulence |
| OprI_PSEAE | Major outer membrane lipoprotein I | 83 | P11221 | Abundant lipoprotein; immunogenic surface antigen |
| OprD_PSEAE | Basic amino acid porin OprD | 429 | P32722 | Carbapenem uptake channel; loss causes carbapenem resistance |
| PilA_PSEAE | Type IV major pilin protein PilA | 150 | P04739 | Major pilin subunit; mediates adhesion and motility |
| PcrV_PSEAE | Type III secretion needle tip protein PcrV | 294 | G3XD49 | T3SS tip protein; key virulence factor and vaccine target |
| ToxA_PSEAE | Exotoxin A | 638 | P11439 | ADP-ribosyltransferase; inhibits protein synthesis |
| LasB_PSEAE | Elastase LasB | 498 | P14756 | Zinc metalloprotease; degrades host tissues and immune factors |
| PhoP_PSEAE | Two-component response regulator PhoP | 225 | Q9I4F9 | Regulator of antimicrobial peptide resistance and virulence |
| AlgD_PSEAE | GDP-mannose 6-dehydrogenase AlgD | 433 | P11759 | Key enzyme in alginate biosynthesis; mucoid phenotype |
| Azurin_PSEAE | Azurin (blue copper protein) | 128 | P00282 | Redox protein with antitumor and immunomodulatory properties |

### 2. `klebsiella_pneumoniae_antigens.fasta`

Antigen and virulence proteins from **Klebsiella pneumoniae** (肺炎克雷伯菌), a Gram-negative pathogen causing pneumonia, urinary tract infections, and sepsis. Hypervirulent (hvKP) and carbapenem-resistant (CRKP) strains are major clinical threats.

| Protein ID | Description | Length (aa) | UniProt/NCBI Accession | Function |
|---|---|---|---|---|
| OmpK36_KLEPN | Porin OmpK36 | 349 | A0A0K1LCY4 | Major porin; porin loss contributes to antibiotic resistance |
| OmpK35_KLEPN | Porin OmpK35 | 362 | A0A2P1E303 | General porin; downregulation associated with resistance |
| OmpA_KLEPN | Outer membrane protein A | 346 | P24017 | Structural protein; involved in adhesion and immune evasion |
| FimH_KLEPN | Type 1 fimbrial adhesin FimH | 300 | B0LF88 | Adhesin mediating host cell attachment; vaccine candidate |
| MrkD_KLEPN | Type 3 fimbrial adhesin MrkD | 295 | P21648 | Adhesin for biofilm formation on medical devices |
| YbtQ_KLEPN | Yersiniabactin siderophore transporter YbtQ | 710 | A4GZC8 | ABC transporter for iron acquisition; virulence determinant |
| EntB_KLEPN | Enterobactin synthase component B | 286 | A0A486V082 | Siderophore biosynthesis; iron scavenging under host restriction |
| KpnO_KLEPN | Outer membrane porin KpnO (OmpK37) | 415 | O87754 | Porin associated with carbapenem resistance phenotype |
| Peg344_KLEPN | DMT family inner membrane transporter Peg-344 | 300 | QXV89895.1 | Marker for hypervirulent K. pneumoniae (pLVPK plasmid) |
| RmpA_KLEPN | Regulator of mucoid phenotype A | 217 | Q8GM82 | Capsule regulator; hypermucoidy in hvKP strains |

## Data Sources

All protein sequences were retrieved from public databases:

- **UniProt Knowledgebase (UniProtKB)**: https://www.uniprot.org
  - Reviewed (Swiss-Prot) entries were preferred when available
  - TrEMBL entries were used for proteins without reviewed sequences
- **NCBI Protein Database**: https://www.ncbi.nlm.nih.gov/protein
  - Used for Peg-344 (accession QXV89895.1) which was not available in UniProt

Sequences correspond to the PAO1 reference strain for *P. aeruginosa* and representative clinical isolates for *K. pneumoniae*.

## FASTA Format

- Header format: `>ProteinID_Species Species_name Protein_description`
- Each sequence is on a single line (no line wrapping)
- Sequences contain only standard amino acid characters (ACDEFGHIKLMNPQRSTVWY)

## Usage

These files can be used as input for the BacterialAntigenFinder pipeline:

```bash
python main.py --input example_data/pseudomonas_aeruginosa_antigens.fasta
python main.py --input example_data/klebsiella_pneumoniae_antigens.fasta
```

## Notes

- Sequences include signal peptides and propeptides where present in the original database entries
- Some sequences (e.g., RmpA) may represent partial-length proteins as full-length sequences were not available in public databases
- These sequences are provided for testing purposes; researchers should verify sequences against current database versions for production use
