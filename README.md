# Microbiome_browser
A shiny app to view and explore microbiome data. 
<p>
This app is designed to demonstrate some basic microbiome analysis techniques on a real-world dataset.The data used in this example was taken from PRJNA822685, a pubically avaliable dataset of oral sqaumous cell carcinoma tumor (abnormal) and normal tumor-adjacent (control) tissue. This publication is avaliable at: https://pubmed.ncbi.nlm.nih.gov/38589501/ and the citation is below:
</p>

<p>
Cai L, Zhu H, Mou Q, Wong PY, Lan L, Ng CWK, Lei P, Cheung MK, Wang D, Wong EWY, Lau EHL, Yeung ZWC, Lai R, Meehan K, Fung S, Chan KCA, Lui VWY, Cheng ASL, Yu J, Chan PKS, Chan JYK, Chen Z. Integrative analysis reveals associations between oral microbiota dysbiosis and host genetic and epigenetic aberrations in oral cavity squamous cell carcinoma.", 
<i>NPJ Biofilms Microbiomes.</i> 2024 Apr 8;10(1):39. doi: 10.1038/s41522-024-00511-x. PMID: 38589501; PMCID: PMC11001959.
</p>

<p>
The dataset was curated to include only tumor and normal tumor tissue from patients with matched samples.
</p>

<p>
To deploy this application, clone the github repo. Use <b>process_reads.Rmd</b> to process the reads (although that data is already avaliable in the file <b>00data.Rds</b>), and run the <b>app.R</b> file from RStudio.
</p>
