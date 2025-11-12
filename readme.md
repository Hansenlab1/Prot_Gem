Proteomics Gem - Analysis Pipeline (PAP)This Shiny application is a complete, end-to-end tool for the analysis of proteomics data. It guides the user through data loading, normalization, and several common downstream analysis methods.The app is built on a modular framework, with each analysis step contained in its own tab and R file.Table of ContentsFeaturesInstallationHow to Run the AppApp Workflow: A Step-by-Step GuideTab 1: Load DataTab 2: NormalizationTab 3: Volcano PlotTab 4: Correlation PlotTab 5: pi0 EstimationTab 6: Marker SignaturesRepository StructureFeaturesModular Pipeline: A 6-tab interface that walks the user from raw data to final analysis.Robust Data Loading: Supports loading separate protein matrix and sample metadata files, ensuring compatibility.Normalization Suite: Includes common normalization methods (Log2, Total Sum, Quantile, Median, etc.) with boxplot visualization.Differential Expression: Generates interactive volcano plots to find up/down-regulated proteins.Correlation Analysis: Correlates protein expression with other proteins or with sample metadata.Statistical Validation: Includes a pi0 estimation tab to assess the proportion of true null hypotheses from the p-value distribution.Advanced Signature Analysis: A powerful module (Tab 6) for complex marker analysis:Marker File Mode: Define complex sum all or ratio signatures in an external file.Dynamic Mode: Explore data by dynamically creating signatures from any annotation column (e.g., "Division", "Category").Meta-Matrix: Run and visualize dozens of comparisons at once in a summary heatmap.Data Source-Switching: Validate all signature analyses by instantly switching between Raw and Normalized data.InstallationTo run this app on your local machine, you will need R and RStudio.Step 1: Install R and RStudioR: Download and install the R programming language from the Comprehensive R Archive Network (CRAN).RStudio: Download and install the free RStudio Desktop IDE from the RStudio (Posit) website.Step 2: Get the App CodeYou can either download this repository as a ZIP file or clone it using Git:git clone \[URL\_TO\_YOUR\_REPO]
Unzip the file if necessary and place the folder in a convenient location.Step 3: Install Required R PackagesOpen RStudio and run the following command in the Console to install all the necessary packages. This may take a few minutes.# Install packages from CRAN
install.packages(c(
"shiny",
"DT",
"ggplot2",
"plotly",
"pheatmap",
"shinycssloaders",
"reshape2",
"dplyr",
"ggrepel"
))

# Install packages from Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install(c("limma", "vsn", "cp4p"))
How to Run the AppOpen RStudio.Go to File > Open Project... and navigate to the folder where you downloaded this repository. Select the .Rproj file (if one exists) or simply open the app.R file.With the app.R file open in the RStudio editor, click the "Run App" button (usually a green arrow) in the top-right corner of the editor pane.The application will launch in a new window, starting on "Tab 1: Load Data".App Workflow: A Step-by-Step GuideThis app is designed to be used in order, from Tab 1 to Tab 6.Tab 1: Load DataThis tab loads all your data. You must load both files to proceed.Logic: This module asks for two separate files. It finds the common samples between the protein data columns and the metadata rows, and validates them. This validated data (raw protein data + sample info) is then passed to all other modules.How to Run (using example data):Click "Browse" for "Upload Protein Matrix" and select Gemcombined\_prot\_insoluble4.csv from the data/ folder.Click "Browse" for "Upload Sample Info" and select sample\_info.csv from the data/ folder.Click the "Load \& Validate Data" button.A validation status and a preview of your sample info will appear. You can now move to Tab 2.Tab 2: NormalizationThis tab applies mathematical transformations to your raw data to correct for technical variation.Logic: The raw protein data from Tab 1 is taken as input. You select a method (e.g., "Quantile", "Total Sum Scaling"). The chosen algorithm is applied, and the resulting normalized data matrix is passed to all downstream analysis tabs. The "None" option simply performs a log2(x+1) transform.How to Run:Select a method from the list.View the boxplots to see the data before and after normalization.The normalized data is now available for Tabs 3-6.Tab 3: Volcano PlotThis tab performs a standard differential expression analysis to find proteins that are significantly different between two groups.Logic: This module uses the normalized data from Tab 2. It identifies all available metadata groups (e.g., "Patient", "Control"). When you select two groups and click "Run", it performs a t-test for every single protein and calculates the p-value and Log2 Fold Change.How to Run:Select "Group 1 (Control)" (e.g., Control).Select "Group 2 (Treatment)" (e.g., Patient).Click "Run Analysis".The plot and tables will populate. The p-values from this analysis are saved for Tab 5.Tab 4: Correlation PlotThis tab lets you explore relationships between variables.Logic: You can perform two types of correlation:Protein vs. Metadata: Correlates the expression of one protein (e.g., "COL1A1") against a numeric metadata variable (e.g., "Age", "Timepoint").Protein vs. Protein: Correlates two different proteins against each other.How to Run:Select the data source (raw or normalized).Choose your correlation type and select the proteins/variables.Click "Generate Plot".Tab 5: pi0 EstimationThis is a statistical validation tab that helps you understand your p-value distribution from Tab 3.Logic: pi0 is the estimated proportion of "true null" hypotheses (i.e., proteins that are not differentially expressed). A high pi0 (e.g., 0.9) is normal. A very low pi0 (e.g., 0.2) might suggest issues with your model or widespread changes. This module takes the list of p-values from the Volcano Plot analysis and uses the cp4p package to estimate pi0.How to Run:First, you must run an analysis in Tab 3.Select an estimation method.Click "Calculate pi0".Tab 6: Marker SignaturesThis is the most advanced analysis tab. It analyzes groups of proteins ("signatures") instead of individual ones.Logic: This module has three separate analysis modes, all of which can be run on Raw or Normalized data (using the "Data Source for Analysis" switch).Single Comparison: Uses an uploaded "Markers File" to define signatures (e.g., total ECM = sum(all Core matrisome proteins)) and compares them between two groups.Meta-Matrix: Runs the "Single Comparison" for many groups at once and shows the results as a heatmap.Dynamic Category Analysis: An exploratory tool. It ignores the "Markers File" and builds signatures dynamically from an "Annotator File" (e.g., it will find all proteins in Division: RBC and make a signature, find all in Category: Collagens and make another).How to Run (Single Comparison):Upload Markers\_module2.csv and StromalAnnotator2.csv from the data/ folder.Select the "Normalized Data" source.Select "Group", "Control", and "Patient".Click "Run Single Analysis".Outputs: This module has its own set of output tabs.Signature List: The main results for the signatures.Full Results Table: Protein-level data for the proteins inside your signatures.Volcano Plot: An interactive plot of the signatures.Heatmap: A heatmap of all proteins in the significant signatures.Feature Viewer: A plot to see the values for a single Gene or Signature.Meta-Matrix: The heatmap output for that specific analysis.Repository StructureProteomics-Analysis-Pipeline/
|
|-- app.R                 (Main app file. Click "Run App" here.)
|-- loader\_module.R       (Code for Tab 1)
|-- normalization\_module.R(Code for Tab 2)
|-- volcano\_module.R      (Code for Tab 3)
|-- correlation\_module.R  (Code for Tab 4)
|-- pi0\_module.R          (Code for Tab 5)
|-- marker\_module.R       (Code for Tab 6)
|-- README.md             (This file.)
|
|-- data/                 (Folder for all example data files.)
|   |-- Gemcombined\_prot\_insoluble4.csv  (Example protein data)
|   |-- sample\_info.csv                (Example sample metadata)
|   |-- Markers\_module2.csv            (Example markers file for Tab 6)
|   |-- StromalAnnotator2.csv          (Example annotator file for Tab 6)
|
|-- images/               (Folder for your screenshots.)
|-- main\_app\_screenshot.png
|-- tab\_1\_loader\_output.png
|-- tab\_2\_normalization\_output.png
|-- tab\_3\_volcano\_output.png
|-- tab\_4\_correlation\_output.png
|-- tab\_5\_pi0\_output.png
|-- tab\_6\_markers\_output.png

