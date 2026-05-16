# VariantXplorer

VariantXplorer is a GUI-based bioinformatics application designed for comprehensive Next-Generation Sequencing (NGS) variant analysis through an interactive and user-friendly interface. The platform integrates multiple stages of genomic analysis, including quality assessment, read trimming, reference mapping, variant calling, filtering, visualization, and automated report generation.

The application is developed using Python and Streamlit, providing an accessible workflow for researchers, students, and bioinformatics users with minimal command-line experience.

---

## Features

- FASTQ quality assessment
- Read trimming and preprocessing
- Reference genome mapping
- Variant calling
- Variant filtering
- Genome and variant visualization
- Automated report generation
- Interactive graphical user interface
- Customizable analysis parameters
- Progress tracking using status/progression bars

---

## Supported Organisms

VariantXplorer supports genomic analysis for multiple organism types, including:

- Human
- Bacteria
- Fungi
- Yeast
- Other organisms with compatible reference genomes and annotation files

---

## Supported Input Files

### Query/Input Files
- FASTQ

### Reference & Annotation Files
- FASTA
- GFF
- GTF
- GFF3

---

## Sequencing Read Support

The application supports both:

- Single-End sequencing reads
- Paired-End sequencing reads

By default, the read type is set to **Single-End**.

Users can switch between read types directly within the graphical interface.

---

## Parameter Customization

Users can modify analysis parameters for each processing stage according to their experimental requirements.

To save modified parameter values:

1. Enter the desired parameter value
2. Press **Enter**
3. Click **OK**

If parameters are not modified, the workflow will automatically execute using the default parameter settings.

---

## Running the Analysis

After:

- Uploading the required files
- Selecting sequencing read type
- Adjusting parameters (optional)

Users can start the workflow by clicking the **Run Analysis** button.

The analysis progress can be monitored using the integrated progression/status bar.

---

## Application Download

The packaged desktop application is available in the **Releases** section of this repository.

Please download:

```text
VariantXplorer.zip
