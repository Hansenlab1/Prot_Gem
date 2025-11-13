# 💎 Proteomics Gem – Analysis Pipeline  
For finding the  gems in your data!
*A modular Shiny application for end-to-end proteomics analysis*

## 📒 Table of Contents
1. [Features](#features)  
2. [Installation](#installation)  
3. [How to Run the App](#how-to-run-the-app)  
4. [App Workflow (Step-by-Step)](#app-workflow-step-by-step)  
   - [Tab 1 – Load Data](#tab-1-load-data)  
   - [Tab 2 – Normalization](#tab-2-normalization)  
   - [Tab 3 – Volcano Plot](#tab-3-volcano-plot)  
   - [Tab 4 – Correlation Plot](#tab-4-correlation-plot)  
   - [Tab 5 – pi0 Estimation](#tab-5-pi0-estimation)  
   - [Tab 6 – Marker Signatures](#tab-6-marker-signatures)  
5. [Repository Structure](#repository-structure)

---

## ✨ Features
- **Modular 6-tab pipeline** guiding users from raw input to high-level signature analysis
- **Robust data ingestion** for protein matrices + sample metadata  
- **Normalization suite** (Log2, Total Sum Scaling, Quantile, Median, etc.)  
- **Differential expression** via interactive volcano plots  
- **Correlation explorer** for protein–protein and protein–metadata relationships  
- **pi0 estimation** for statistical validation of p-value distributions  
- **Advanced signature analysis**, including dynamic, marker-file, and meta-matrix modes  
- **Raw vs Normalized data toggling** for validation and comparison  

---

## 💻 Installation

### 1. Install R and RStudio
Download and install:
- **R** from CRAN  
- **RStudio** from Posit  

### 2. Download or Clone the Repository
```bash
git clone https://github.com/yourname/Proteomics-Analysis-Pipeline.git
```

### 3. Install Required Packages
```r
install.packages(c(
  "shiny", "DT", "ggplot2", "plotly", "pheatmap",
  "shinycssloaders", "reshape2", "dplyr", "ggrepel"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma", "vsn", "cp4p"))
```

---

## 🚀 How to Run the App
1. Open **RStudio**  
2. Open the project folder  
3. Open `app.R`  
4. Click **Run App**

The app opens in a new window starting on **Tab 1: Load Data**.

or..

type in the R console: 
> shiny::runApp(launch.browser = TRUE)

---

# 🧭 App Workflow (Step-by-Step)

---

## 🔹 Tab 1 – Load Data

This module:
- Accepts a **Protein Matrix** CSV file (metadata rows need a header with "_meta")  
- Validates format and sample name overlap  
- Passes validated objects to all downstream modules  

### Example Steps:
1. Upload `Gemcombined_prot_insoluble4.csv`  
2. Click **Load & Validate Data**

### Screenshot Placeholder
![Tab 1 Loader Output](images/tab_1_loader_output.png)

---

## 🔹 Tab 2 – Normalization

This module:
- Performs normalization of protein data  
- Displays before/after boxplots  
- Makes normalized matrix available to Tabs 3–6  

Supported methods:
- Log2 transform (default)  
- Quantile  
- Total Sum Scaling  
- Median  
- None (log2(x+1) only)
- SVN
- Loess

### Screenshot Placeholder
![Tab 2 Normalization Output](images/tab_2_normalization_output.png)

---

## 🔹 Tab 3 – Volcano Plot

This module runs:
- Group-vs-group t-tests  
- Computes log2 fold changes  
- Generates interactive volcano plots  
- Saves the p-values for usage in Tab 5  

### How to Run:
1. Select **Group 1**  
2. Select **Group 2**  
3. Click **Run Analysis**
4. you can download the results, or just significantly up- and down- results for further Gene Ontology analysis

### Screenshot Placeholder
![Tab 3 Volcano Output](images/tab_3_volcano_output.png)

---

## 🔹 Tab 4 – Correlation Plot

This module supports:
- **Protein vs Protein** correlations (not yet!) 
- **Protein vs Metadata** correlations  

Data sources:
- Raw  
- Normalized  

### How to Run:
1. Select **Raw** or **Normalized**  
2. Choose correlation type  
3. Select variables  
4. Click **Generate Plot**
5. you can "download"" a svg or pdf version of the plot
6. "download" significant positive and negative correlations for further data analysis

### Screenshot Placeholder
![Tab 4 Correlation Output](images/tab_4_correlation_output.png)

---

## 🔹 Tab 5 – π0 (pi0) Estimation

This module:
- Uses p-values from Tab 3  
- Estimates the proportion of true null hypotheses  
- Helps validate whether observed differences are global vs targeted  

Estimation methods include: (these need testing!)
- Bootstrap  
- Smoother  
- Histogram-based approaches  

### Screenshot Placeholder
![Tab 5 pi0 Output](images/tab_5_pi0_output.png)

---

## 🔹 Tab 6 – Marker Signatures

### Three Analysis Modes

#### **1. Single Comparison Mode**
- Requires **Markers File** + **Annotator File**  
- Computes multi-protein signatures (sum, ratio, composite signatures)  
- Compares between two metadata groups  

#### **2. Meta-Matrix Mode**
- Runs multiple comparisons at once  
- Outputs a heatmap summarizing comparisons  

#### **3. Dynamic Category Analysis**
- Ignores marker files  
- Automatically constructs signatures from annotation categories  
  (e.g., all “Collagens,” all “Immune Proteins,” etc.)

### How to Run (Single Comparison):
1. Upload `Markers_module2.csv`  
2. Upload `StromalAnnotator2.csv`  
3. Select analysis mode  
4. Choose groups  
5. Click **Run Single Analysis**

### Screenshot Placeholder
![Tab 6 Marker Output](images/tab_6_markers_output.png)

---

# 📂 Repository Structure

```
Proteomics-Analysis-Pipeline/
│
├── app.R
├── loader_module.R
├── normalization_module.R
├── volcano_module.R
├── correlation_module.R
├── pi0_module.R
├── marker_module.R
├── README.md
│
├── data/
│   ├── Gemcombined_prot_insoluble4.csv
│   ├── Markers_module2.csv
│   ├── StromalAnnotator2.csv
│
└── images/
    ├── main_app_screenshot.png
    ├── tab_1_loader_output.png
    ├── tab_2_normalization_output.png
    ├── tab_3_volcano_output.png
    ├── tab_4_correlation_output.png
    ├── tab_5_pi0_output.png
    ├── tab_6_markers_output.png
```

---