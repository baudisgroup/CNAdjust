# CNAdjust

## Prerequisites

1. Java 11 or later
2. Local install (recommended) or [docker](https://www.docker.com/)  

For locally installing all dependencies, [R](https://www.r-project.org/) is required and then run with the following R code:

```R
install.packages(c('R.utils','MASS','remotes'))
remotes::install_github('baudisgroup/labelSeg')
```

## Installation

1. Install nextflow by using the follwowing command:

```bash
curl -s https://get.nextflow.io | bash
```

This command will create a `nextflow` file in the current directory.

2. Make Nextflow executable:

```bash
chmod +x nextflow
```

3. Install the pipeline

The first option is to use nextflow, it will be installed in the `$HOME` directory under the `.nextflow/assets` sub-directory.

```bash
./nextflow pull hangjiaz/CNAdjust
```

The second option is to clone the pipeline repository into a desired directory. 

```bash
git clone https://github.com/hangjiaz/CNAdjust.git
```

## Usage

```bash
./nextflow run hangjiaz/CNAdjust -r main --inputdir /path/to/Inputdir --series <series1>,<series2> --outputdir /path/to/Outputdir
```

## Input files

The input files are organized by series, with each series containing the corresponding segment data and associated external information. Each series forms one channel in Nextflow and is processed in parallel.

```
Inputdir/
├── <series1>
│   ├── <series1>.seg.txt            segment data with callings
│   ├── cohort-assignment.txt        cohort assignments
│   ├── cohort-cna-pattern.txt       cohort-specific CNA occurrence
│   ├── sampleid-mapping.txt         optional: sample ID mapping
│   └── cohort-cna-region.txt        optional: custom genomic bin location
└── <series2>
    └── ...     
```

### 1. Segment files (required)

These files should be tab-delimited and contain at least six columns, including sample ID, chromosome, start, end, number of markers, and log R ratio. Additionally, a column named "label" should be included to store CNA callings as strings ("+1" for gain, "-1" for loss, "0" for copy neutral). For level-specific callings, you can use "+2", "+1", "0", "-1", "-2" to represent high-level duplication, low-level duplication, copy neutral, low-level deletion, and high-level deletion respectively. To be recognized by the workflow, segment file names should include ".seg" text.

### External information

External information includes cohort assignments for samples and cohort-specific CNA prior patterns. 

If external information is not provided, you should set the parameter `use_external_ref` to "false". This setting assumes that all samples within the series are biologically and technologically similar. Consequently, the workflow utilizes the CNA frequency calculated from all samples as prior information.

### 2. Cohort assignment

The cohort assignment file should be tab-delimited and include two columns: sample IDs and cohort identifiers. The default filename for this file is "cohort-assignment.txt" and can be changed by the parameter `cohort_assign_file`.

### 3. Cohort CNA occurence 

The cohort CNA occurrence file should also be tab-delimited. Different columns represent CNA occurence probability from different cohorts. Column names should be consistent with cohort identifiers specified in the cohort assignment file. More details about CNA occurence probability see [below](#Cohort-CNA-pattern-data-preparation). The default filename for this file is "cohort-cna-pattern.txt" and can be changed by the parameter `prior_file`.

### 4. Sample ID mapping file (optional)

This optional file is used for mapping IDs used in segment data to expected sample identifiers, such as mapping UUIDs to barcodes in TCGA data. If provided and the parameter `use_idmapping` is set to "true", the sample IDs used in the cohort assignment file should be the expected identifiers. All sample identifiers used in output files will be mapped to these new identifiers. The file format should be tab-delimited and include two columns: original sample identifiers and new sample identifiers. The default filename for this file is "sampleid-mapping.txt" and can be changed by the parameter `idmapping_file`.

### 5. Genomic region file (optional)

This optional file indicates the genomic regions used to calculate cohort CNA occurrence. If provided and the parameter `use_custom_region` is set to "true", the prior computation will be based on the provided regions. The file format should be tab-delimited and include four columns: region index, chromosome, start position, and end position (example see `data/hg38_bin.txt`). The default filename for this file is "cohort-cna-region.txt" and can be changed by the parameter `region_file`.

## Output

By default, the output directory is a folder named "output" located in the current directory. The structure of the output files is as follows:

```
Outputdir/
├── <series1>
│   ├─── data_quality_report.txt
│   ├─── <sample1>
│   │   ├── segments.pdf
│   │   └── result.seg.txt
│   ├── <sample2>
│   │   ├── segments_before_shift.pdf
│   │   ├── segments.pdf
│   │   └── result.seg.txt
│   └── ...
└─ <series2>   
    └── ...                 
```

1. **Series Folder**: Each series folder contains a data quality report summarizing the features extracted from segment profiles and the manipulations performed on them.

    Manipulations on the data include:

    * "initial": No change to the data.
   
    * "shift-baseline-lower/higher": Attempt to shift the baseline, with the result being better than the original one, therefore used as the final output.
   
    * "shift-baseline-lower/higher_reverse": Attempt to shift the baseline, but the result was not satisfactory, so the original callings were retained.
   
    * "noisy_lowthre_\<lowcutoff\>\_highthre\_\<highcutoff\>": Using cutoffs for noisy samples.
   
    * "noisy_nochange": Attempt to apply cutoffs for noisy samples, but the result was not satisfactory, so the original callings were retained.

3. **Sample Folders**: Within each series folder, there are individual folders for each sample. These folders contain:

    Segment Data: Final segment data files with adjusted callings. The default filename for this file is "result.seg.txt" and can be changed by the parameter `output_seg`.
   
    Segment Plots: Visual representations of segments colored by CNA states. The default filename for this file is "segments.pdf" and can be changed by the parameter `output_plot`. If baseline shifting was performed for the sample, an additional plot named "<plotfileBasename>_before_shift.pdf" is included to show the original callings before shifting.

    
## Parameters

Mandatory parameters

* `--inputdir`: Path to the input data directory.
* `--series`: Name of the series to be analyzed. It should be the folder name of the subfolder under `inputdir`.

Optional parameters

* `--outputdir`: Path to the output data directory.
* `--use_external_ref`: Logical value to determine if external information is used to adjust CNA.
* `--use_idmapping`: Logical value to determine if sample identidier mapping is applied  given the sampleid mapping file.
* `--use_custom_region`: Logical value to determine if a custom prior is applied given the genomic bin location file.
* `--genome`: Reference genome assembly version used in segment data. Available options are "hg38", "hg19", and "hg18".
* `--cohort_assign_file`: Filename of the input cohort assignment file.
* `--prior_file`:Filename of the input cohort CNA occurrence file.
* `--idmapping_file`: Filename of the input sample id mapping file.
* `--region_file`: Filename of the input genomic region file.
* `--output_seg`: Filename of the output segment data.
* `--output_plot`: Filename of the output plot of segment data.
* `--logrsd`: Upper limit of weighted standard deviation of logR by the number of involved markers to identify potential problematic samples in the strict criteria.
* `--cnafrac`: Upper limit of CNA fraction to identify potential problematic samples in the strict criteria.
* `--dupdelratio`: Upper limit of the fraction ratio and its reciprocal as lower limit to identify potential problematic samples in the strict criteria.
* `--cnafrac2`: Upper limit of CNA fraction to identify potential problematic samples in the more relaxed criteria 2.
* `--dupdelratio2`: Upper limit of the fraction ratio and its reciprocal as lower limit to identify potential problematic samples in the more relaxed criteria 2.
* `--cnafrac3`: Upper limit of CNA fraction to identify potential problematic samples in the more relaxed criteria 3.
* `--segnum`: Upper limit of segment number to identify noisy profiles.
* `--lowthre`: Thresholds to call low-level CNA in noisy profiles. It should be a positive value used to call duplication, and its opposite is used to call deletion.
* `--highthre`: Thresholds to call high-level CNA in noisy profiles. It should be a positive value used to call duplication, and its opposite is used to call deletion.

## Profile 

By default, the pipeline is locally executed, but it can be run using docker engine by simply setting `-profile docker`.

## Cohort CNA pattern data preparation

By default, the prior CNA pattern in each cohort represents the occurrence probability of CNAs in genomic bins. The bin size is 1 MB. For GRCh38, there are 2 * 3106 bins (gain_frequency + loss_frequency). Detailed genomic bin locations for GRCh38, GRCh37, and GRCh36 can be found in `data/hg38_bin.txt`, `data/hg19_bin.txt`, and `data/hg18_bin.txt` respectively. You can define the CNA prior pattern yourself based on the same genomic region or custom genomic regions by setting the parameter `--use_custom_region` to "true", as long as the values can be used as probabilities. Here we introduce two simple approaches to get the cohort-specific CNA reference pattern.

### from Progenetix

The simplest way to get the cohort CNA pattern is to access pre-computed CNA frequency data from the Progenetix database via the Beacon v2 API. It is defined as the percentage of samples showing a CNV for a genomic region over the total number of samples in a cohort specified by filters. For example, if you want to get the reference CNA pattern from glioblastoma samples, you can use the pgxRpi R package to query:

```r
library(pgxRpi)
library(SummarizedExperiment)

# choose the filter to specify glioblastoma cohort using NCIt encoding system
tumor_filter <- "NCIT:C3058" 
# query
cnafreq <- pgxLoader(type = 'frequency',filters = tumor_filter,output = 'pgxmatrix')  
# transform from frequency value to probability value
cnaprob <- assay(cnafreq)/100 
# Write data in the format that the workflow accepts
write.table(cnaprob,"cohort-cna-pattern.txt"),quote = F,sep = '\t',row.names = F,col.names = T)
```

You can also get the CNA frequency data directly via the [REST API](https://docs.progenetix.org/file-formats/?h=#data-matrix-files), but the data need to be modified slightly to follow the format which CNAadjust uses.

### from custom segment data 

You can utilize the function `segtoFreq` from `pgxRpi` R package (version >= 1.0.1 or >= 1.1.2) to derive CNA frequency from segment data. These frequencies need to be transformed into probabilities. By default, the binning aligns with that used in CNAdjust. However, if you calculate CNA frequency for different genomic bins, such as using different bin sizes, you should provide the [region file](#5-Genomic-region-file-optional) accordingly. Example usage is as follows:

```bash
nextflow run hangjiaz/CNAdjust --inputdir /path/to/Inputdir --series <series1>,<series2> --outputdir /path/to/Outputdir --use_custom_prior true 
```

