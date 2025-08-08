<div align="center">
   
  <h1>Microscope🔬</h1>
  
  **_Microbiome Diversity and Taxonomy Explorer_**

  [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/microscope?color=brightgreen)](https://github.com/OpenOmics/microscope/issues) [![GitHub license](https://img.shields.io/github/license/OpenOmics/microscope)](https://github.com/OpenOmics/microscope/blob/main/LICENSE) 
  
  <i>
    This is the home of the shiny application, microscope. Its long-term goals are to empowers researchers to rapidly explore, compare, and visialize microbiome diversity and taxonomy.
  </i>
</div>


## Overview

`microscope` is an [R Shiny applications](https://shiny.posit.co) that provides an interactive interface for visualizing **diversity metrics** and **taxonomy abundance** from microbiome data. Users can upload metadata, diversity results and taxonomy abundance files, select grouping variables, and generate publication-ready plots.

Microscope’s mission is to accelerate microbiome insights by providing a user-friendly interface for exploring and visualizing microbiome data. It aims to bridge the gap between raw microbiome outputs and interpretable results, enabling researchers to rapidly explore, compare, and visualize microbiome diversity and taxonomy.

### Features

Here is a list of features that microscope provides:

- Upload metadata (`.csv`) and specify the folder containing diversity files (compatible with metamorph metamorph-generated diversity metrics).
- Dynamically select grouping variables for visualization.
- Generate **violin plots**, **boxplots**, and **scatter plots** for alpha diversity metrics.
- Generate **PCoA plots** for beta diversity metrics.
- Optional statistical comparisons between groups.
- Taxonomy stacked bar plots

## Getting Started

### Dependencies

**Requires:** `R>=4.0`

Before getting started, please ensure you have the following dependencies installed:

- [R](https://cloud.r-project.org): To install any required R packages and run the Shiny app.

### Installation

Once you have installed the dependencies above, please run the following commands to setup and install `microscope`:

```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/microscope.git
# Change your working directory
cd microscope/
# Assume you have R installed, 
# install pak to install the
# any required R packages
Rscript -e "if (!requireNamespace('pak', quietly = TRUE)) install.packages('pak', repos = 'https://cloud.r-project.org')"
# Install the required R packages
# from the DESCRIPTION file
Rscript -e "pak::local_install_dev_deps('.', dependencies = TRUE)"
```

### Running microscope

Once you have installed the required packages, you can run the shiny application using the following command:

```bash
# Run the shiny app on your local machine
Rscript -e "shiny::runApp('shiny/', launch.browser = TRUE)"
```