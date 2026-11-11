# **Data Preprocessing**

## **1\. Overview**

Preprocessing in PHANES is designed to handle:

* Heterogeneous biological features  
* Length-dependent biases  
* Missing experimental annotations

Such preprocessing pipelines are standard in bioinformatics and proteomics analyses, where normalization, transformation, and imputation are required for meaningful comparisons across proteins [(1–4)](https://paperpile.com/c/ImpoB4/8QcL+jpq0+8dgO+9WtY).

---

## **2\. Length Normalization in Protein Data**

Protein features are often length-dependent, particularly aggregated properties.

All features are checked for association with length via Linear Regression modeling and 

### 2a. Normalization by sequence length is applied to extensive features:

[![x\_{norm} = x / L][image1]](https://saxarona.github.io/mathjax-viewer/?input=x_%7Bnorm%7D%20%3D%20x%20%2F%20L#0) [(5)](https://paperpile.com/c/ImpoB4/g4sm)

Where:

* [![x = feature\\, value][image2]](https://saxarona.github.io/mathjax-viewer/?input=x%20%3D%20feature%5C%2C%20value#0)  
* [![L = protein\\, length][image3]](https://saxarona.github.io/mathjax-viewer/?input=L%20%3D%20protein%5C%2C%20length#0)

### 2b. Normalization by Residuals 

Features that power-law dependence to length are normalized using allometric residualization. Log-log regressions of specific features against protein length were fitted and the residuals were used as size-corrected features.[(6–8)](https://paperpile.com/c/ImpoB4/rQC4+oFlK+8ipa)

[![Residuals = \\log(Y) - \\log(Y\_{pred})][image4]](https://saxarona.github.io/mathjax-viewer/?input=Residuals%20%3D%20%5Clog\(Y\)%20-%20%5Clog\(Y_%7Bpred%7D\)#0)

---

## **3\. Transformation of Biological Features**

Biological datasets frequently exhibit:

* Skewness  
* Heavy-tailed distributions

The Yeo–Johnson transformation is applied to stabilize variance [(9,10)](https://paperpile.com/c/ImpoB4/13q1+aQlW):

For x ≥ 0:  
[![y = ((x + 1)^{\\lambda} -1) / \\lambda][image5]](https://saxarona.github.io/mathjax-viewer/?input=y%20%3D%20\(\(x%20%2B%201\)%5E%7B%5Clambda%7D%20-1\)%20%2F%20%5Clambda#0)

For x \< 0:  
[![y = -\[((-x + 1)^{2 - \\lambda} -1) / (2 - \\lambda)\]][image6]](https://saxarona.github.io/mathjax-viewer/?input=y%20%3D%20-%5B\(\(-x%20%2B%201\)%5E%7B2%20-%20%5Clambda%7D%20-1\)%20%2F%20\(2%20-%20%5Clambda\)%5D#0)

---

## **4\. Standardization (Z-score Normalization)**

Features are standardized using [(11)](https://paperpile.com/c/ImpoB4/oTP2):

[![z = (x - \\mu) / \\sigma][image7]](https://saxarona.github.io/mathjax-viewer/?input=z%20%3D%20\(x%20-%20%5Cmu\)%20%2F%20%5Csigma#0)

Where:

* [![\\mu = mean][image8]](https://saxarona.github.io/mathjax-viewer/?input=%5Cmu%20%3D%20mean#0)  
* [![\\sigma = standard\\, deviation][image9]](https://saxarona.github.io/mathjax-viewer/?input=%5Csigma%20%3D%20standard%5C%2C%20deviation#0)

This ensures comparability across features.

---

## **5\. Missing Data Imputation**

Missing values (in thermostability and abundance features) are handled using KNN imputation[(12,13)](https://paperpile.com/c/ImpoB4/CoMT+7fD2):

[![x\_{missing} = average\\, of \\, k\\, nearest\\, neighbors\\, (k = 15)][image10]](https://saxarona.github.io/mathjax-viewer/?input=x_%7Bmissing%7D%20%3D%20average%5C%2C%20of%20%5C%2C%20k%5C%2C%20nearest%5C%2C%20neighbors%5C%2C%20\(k%20%3D%2015\)#0)

This preserves local feature structure.

---

## **6\. Min-Max Scaling (Normalization)**

Per residue features are scaled with Min-Max scaling using [(11)](https://paperpile.com/c/ImpoB4/oTP2):

[![X\_{scaled} = (X - X\_{min}) / (X\_{max} - X\_{min})][image11]](https://saxarona.github.io/mathjax-viewer/?input=X_%7Bscaled%7D%20%3D%20\(X%20-%20X_%7Bmin%7D\)%20%2F%20\(X_%7Bmax%7D%20-%20X_%7Bmin%7D\)#0)

This sets the values in the \[0,1\] range which already is the case for the Physical/Chemical features. In order to be able to compare direct values for localized features we need the range of values to be similar across all per residue features.

---

# **Methods**

## **1\. Dimensionality Reduction Methods**

### **1.1 Principal Component Analysis (PCA)**[(14)](https://paperpile.com/c/ImpoB4/3mv3)

PCA reduces dimensionality by projecting data onto uncorrelated principal components while retaining as much as possible of the variation present in the dataset Preprocessing of data is important because of the dominance of features with the largest variances in the Principal components[(14)](https://paperpile.com/c/ImpoB4/3mv3):

[![Z = XW][image12]](https://saxarona.github.io/mathjax-viewer/?input=Z%20%3D%20XW#0)

Where:

* [![X = standardized\\,feature\\,matrix][image13]](https://saxarona.github.io/mathjax-viewer/?input=X%20%3D%20standardized%5C%2Cfeature%5C%2Cmatrix#0)  
* [![W = eigenvectors\\,of\\,covariance\\,matrix][image14]](https://saxarona.github.io/mathjax-viewer/?input=W%20%3D%20eigenvectors%5C%2Cof%5C%2Ccovariance%5C%2Cmatrix#0)

Data Preprocessing: Length Normalization, Transformation, KNN Imputation, Standarization

---

### **1.2 UMAP (Uniform Manifold Approximation and Projection)**

UMAP is a dimensionality reduction algorithm that captures non-linear relationships in high-dimensional data that can be plotted similarly to PCA[(15,16)](https://paperpile.com/c/ImpoB4/kwNI+udIE). It is ideal as a general purpose dimensionality reduction method and it yields similar results to t-SNE[(17)](https://paperpile.com/c/ImpoB4/uvA7) in single-cell data but it’s significantly faster [(16,18)](https://paperpile.com/c/ImpoB4/udIE+5APh)

Distance metric used [(18)](https://paperpile.com/c/ImpoB4/5APh):

Data Preprocessing: Length Normalization, Transformation, KNN Imputation, Standarization

---

## **2\. Plot-Specific Methods**

### **2.1 Spider Plot (Z-score Representation)**

Spider plots use Z-scores to compare and identify the most divergent protein features against the Human Proteome or a user-defined set and plots in terms of standard deviation. Proteins tend to have similar trends across the same categories of physical/chemical properties, for this reason once such a feature is selected, we eliminate all other features from the same category. In order to further reduce redundancy we exclude features that strongly correlate with other selected features (Spearman’s Correlation Coefficient \> 0.9).

[![z = (x - \\mu) / \\sigma][image15]](https://saxarona.github.io/mathjax-viewer/?input=z%20%3D%20\(x%20-%20%5Cmu\)%20%2F%20%5Csigma#0)

Data Preprocessing: Length Normalization, Transformation, KNN imputation

---

### **2.2 Heatmap Clustering**

Hierarchical clustering (e.g., Ward’s method) is the process of classification of unclassified elements into agnostic groups. It is widely used in taxonomy and retrieval of information even for tasks where a non hierarchical grouping is needed. The grouping used in this application is based on Euclidean distance [(19)](https://paperpile.com/c/ImpoB4/RNVB). In order to portray different scale features in the same heatmap data preprocessing is implemented:

[![d(x,\\,y) = \\sqrt{\\sum (x\_i - y\_i)^2}][image16]](https://saxarona.github.io/mathjax-viewer/?input=d\(x%2C%5C%2Cy\)%20%3D%20%5Csqrt%7B%5Csum%20\(x_i%20-%20y_i\)%5E2%7D#0)

Data Preprocessing: Length Normalization, Transformation, KNN Imputation, Standarization

---

### **2.3 Protein Similarity Network**

Graph networks are quite practical to depict protein similarity networks in many cases [(20–23)](https://paperpile.com/c/ImpoB4/TzmI+xF0Y+ZnMW+9v33). Here a graph network is illustrated, where selected proteins are nodes and edges are included when similarity \> 0.5 for any 2 proteins[(20)](https://paperpile.com/c/ImpoB4/TzmI) .

Pairwise protein similarity is quantified with cosine similarity, that unlike euclidean distance, is invariant to the overall magnitude of a protein’s feature vector and therefore compares the angles in feature space as used widely in scientific literature [(24–26)](https://paperpile.com/c/ImpoB4/zIIm+gQAF+hWVv). Because cosine similarity is sensitive to the input representation and can yield uninformative score when preprocessing is neglected [(27)](https://paperpile.com/c/ImpoB4/0Src), features are handled accordingly before comparison by excluding features with zero or near zero variance [(28)](https://paperpile.com/c/ImpoB4/AGP1), then we follow the length normalization, transformation, KNN imputation, standard scaling as described and omit highly collinear features. Although cosine similarity is not translation invariant, centering removes the large shared offset of the predominantly non-negative descriptors, which would otherwise drive all pairwise similarities towards 1, yielding the angular similarity of each protein’s deviation from the population mean, which is equivalent to the Pearson correlation coefficient [(29)](https://paperpile.com/c/ImpoB4/xDZc)

[![cosine\\,similarity = (x \\cdot y) / (\\parallel x\\parallel\\,\\parallel y \\parallel)][image17]](https://saxarona.github.io/mathjax-viewer/?input=cosine%5C%2Csimilarity%20%3D%20\(x%20%5Ccdot%20y\)%20%2F%20\(%5Cparallel%20x%5Cparallel%5C%2C%5Cparallel%20y%20%5Cparallel\)#0)

Data Preprocessing: Length Normalization, Transformation, KNN Imputation, Standarization

---

### **2.4 Violin Plot (Distribution Estimation)**

Violinplots are informative tools for displaying feature distributions, which are estimated using kernel density estimation (KDE). Similar to the boxplot, contains information about the median, quartiles, IQR and outliers [(30–32)](https://paperpile.com/c/ImpoB4/KiI9+GNeu+K5Hn):

[![f(x) = (1 / nh) \\sum K((x - x\_i)/h)][image18]](https://saxarona.github.io/mathjax-viewer/?input=f\(x\)%20%3D%20\(1%20%2F%20nh\)%20%5Csum%20K\(\(x%20-%20x_i\)%2Fh\)%20#0)

Data Preprocessing: None

---

### **2.5 Protein Group Comparison**

Additionally, PHANES enables protein group comparison by effect size. Effect size explains the impact of a feature in the distinction of two groups in comparison with the use of Cohen’s d parameter [(33)](https://paperpile.com/c/ImpoB4/8xMS). Hedge proved that for small sample sizes, Cohen’s d is biased and proposed a correction parameter based on the degrees of freedom.Durlak later gave a simpler approximation of this factor [(34,35)](https://paperpile.com/c/ImpoB4/80FL+o3aH). These metrics assume normality of distributions. This implementation identifies the 10 features with the highest effect size metric, so in case there is a tie, tiebreaker role is assumed by relative difference. 

#### **Cohen’s d (Effect Size) if n \>=50** [(33,35)](https://paperpile.com/c/ImpoB4/8xMS+o3aH)

[**![d = (\\mu\_1 - \\mu\_2) / s\_{pooled}][image19]**](https://saxarona.github.io/mathjax-viewer/?input=d%20%3D%20\(%5Cmu_1%20-%20%5Cmu_2\)%20%2F%20s_%7Bpooled%7D#0)

[![s\_{pooled} = \\sqrt{((n\_1 - 1)s\_1^2 + (n\_2-1)s\_2^2) / (n\_1 + n\_2 -2)}][image20]](https://saxarona.github.io/mathjax-viewer/?input=s_%7Bpooled%7D%20%3D%20%5Csqrt%7B\(\(n_1%20-%201\)s_1%5E2%20%2B%20\(n_2-1\)s_2%5E2\)%20%2F%20\(n_1%20%2B%20n_2%20-2\)%7D%20#0)

#### **Hedge’s g if n \<50** [(34)](https://paperpile.com/c/ImpoB4/80FL)

Correction parameter(Hedge): [![(1 - 3/(4n-9)), n = n\_1+n\_2][image21]](https://saxarona.github.io/mathjax-viewer/?input=\(1%20-%203%2F\(4n-9\)\)%2C%20n%20%3D%20n_1%2Bn_2#0)  
Correction parameter(Durlak): [![n - 3 / n - 2.25 \* \\sqrt{(n-2)/n}][image22]](https://saxarona.github.io/mathjax-viewer/?input=n%20-%203%20%2F%20n%20-%202.25%20*%20%5Csqrt%7B\(n-2\)%2Fn%7D#0)  
[![g = d \* Correction\\, parameter(Hedge)][image23]](https://saxarona.github.io/mathjax-viewer/?input=g%20%3D%20d%20*%20Correction%5C%2C%20parameter\(Hedge\)#0)

Data Preprocessing: Transformation

#### **Relative Difference**

[**![Relative\\, Difference = \\mid\\mu\_1 - \\mu\_2\\mid / (\\mid\\mu\_1\\mid + \\mid\\mu\_2\\mid) / 2][image24]**](https://saxarona.github.io/mathjax-viewer/?input=Relative%5C%2C%20Difference%20%3D%20%5Cmid%5Cmu_1%20-%20%5Cmu_2%5Cmid%20%2F%20\(%5Cmid%5Cmu_1%5Cmid%20%2B%20%5Cmid%5Cmu_2%5Cmid\)%20%2F%202#0)

---

### **2.6 Profile Plot Scaling**

The latest addition are scatterplots which illustrate per residue features. This way the user is able to directly analyse protein sequences and enables consensus/comparison for specific parts of interest in proteins.

Data Preprocessing: Min-Max Scaling [(11)](https://paperpile.com/c/ImpoB4/oTP2)

---

### **2.7 Feature Redundancy**

In this ensemble we provide pairwise Spearman’s ρ[(36)](https://paperpile.com/c/ImpoB4/g1sV) heatmap that indicates the pairwise monotonic correlation between any 2 features. This serves as an indicative first step of identifying relationships between features and **not as a complete thorough feature redundancy implementation**. For this implementation only the scalar values are used,  **not the per-residue feature values**. The 0.9 threshold mentioned is approximated on the basis of its correspondence to a Variance Inflation Factor (VIF) \> 5 [(37)](https://paperpile.com/c/ImpoB4/KacF) under rank transformation. This threshold is not directly indicative of redundant features because it does not refer to multicollinearity. 

Data Preprocessing: None

---

1\.	[Papadakis S, Garefalakis A, Lemonakis C, Chimonaki C, Zopounidis C. Machine Learning Applications for Accounting Disclosure and Fraud Detection \[Internet\]. IGI Global; 2020\. 291 p. Available from: https://play.google.com/store/books/details?id=19IPEAAAQBAJ](http://paperpile.com/b/ImpoB4/8QcL)

2\.	[Hsu LL, Culhane AC. Impact of Data Preprocessing on Integrative Matrix Factorization of Single Cell Data. Front Oncol \[Internet\]. 2020 June 23;10:973. Available from: http://dx.doi.org/10.3389/fonc.2020.00973](http://paperpile.com/b/ImpoB4/jpq0)

3\.	[Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol \[Internet\]. 2019 Dec 23;20(1):296. Available from: http://dx.doi.org/10.1186/s13059-019-1874-1](http://paperpile.com/b/ImpoB4/8dgO)

4\.	[Bolstad BM, Irizarry RA, Astrand M, Speed TP. A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. Bioinformatics \[Internet\]. 2003 Jan 22;19(2):185–93. Available from: http://dx.doi.org/10.1093/bioinformatics/19.2.185](http://paperpile.com/b/ImpoB4/9WtY)

5\.	[Brocchieri L, Karlin S. Protein length in eukaryotic and prokaryotic proteomes. Nucleic Acids Res \[Internet\]. 2005 June 10;33(10):3390–400. Available from: http://dx.doi.org/10.1093/nar/gki615](http://paperpile.com/b/ImpoB4/g4sm)

6\.	[Adiguzel Y. Information-theoretic approach in allometric scaling relations of DNA and proteins. Chem Biol Drug Des \[Internet\]. 2022 Feb;99(2):331–43. Available from: http://dx.doi.org/10.1111/cbdd.13988](http://paperpile.com/b/ImpoB4/rQC4)

7\.	[Oestereicher MA, da Silva-Buttkus P, Gailus-Durner V, Marschall S, Fuchs H, IMPC Consortium, et al. Rethinking ratio-based normalization towards model-based approaches in heart weight analysis. Sci Rep \[Internet\]. 2026 Mar 17;16(1). Available from: http://dx.doi.org/10.1038/s41598-026-43503-x](http://paperpile.com/b/ImpoB4/oFlK)

8\.	[Liu L, Hyeon C. Contact Statistics Highlight Distinct Organizing Principles of Proteins and RNA. Biophys J \[Internet\]. 2016 June 7;110(11):2320–7. Available from: http://dx.doi.org/10.1016/j.bpj.2016.04.020](http://paperpile.com/b/ImpoB4/8ipa)

9\.	[Website \[Internet\]. Available from: https://doi.org/10.1093/biomet/87.4.954](http://paperpile.com/b/ImpoB4/13q1)

10\.	[Zwiener I, Frisch B, Binder H. Transforming RNA-Seq data to improve the performance of prognostic gene signatures. PLoS One \[Internet\]. 2014 Jan 8;9(1):e85150. Available from: http://dx.doi.org/10.1371/journal.pone.0085150](http://paperpile.com/b/ImpoB4/aQlW)

11\.	[Bishop CM. Pattern Recognition and Machine Learning \[Internet\]. Springer New York; \[cited 2026 May 4\]. 20 p. Available from: https://link.springer.com/book/9780387310732](http://paperpile.com/b/ImpoB4/oTP2)

12\.	[Di Lena P, Sala C, Prodi A, Nardini C. Missing value estimation methods for DNA methylation data. Bioinformatics \[Internet\]. 2019 Oct 1;35(19):3786–93. Available from: http://dx.doi.org/10.1093/bioinformatics/btz134](http://paperpile.com/b/ImpoB4/CoMT)

13\.	[Lazar C, Gatto L, Ferro M, Bruley C, Burger T. Accounting for the Multiple Natures of Missing Values in Label-Free Quantitative Proteomics Data Sets to Compare Imputation Strategies. J Proteome Res \[Internet\]. 2016 Apr 1;15(4):1116–25. Available from: http://dx.doi.org/10.1021/acs.jproteome.5b00981](http://paperpile.com/b/ImpoB4/7fD2)

14\.	[Jolliffe IT. Principal Component Analysis \[Internet\]. Springer New York; \[cited 2026 May 3\]. 30 p. Available from: https://link.springer.com/book/10.1007/b98835](http://paperpile.com/b/ImpoB4/3mv3)

15\.	[McInnes L, Healy J, Saul N, Großberger L. UMAP: Uniform Manifold Approximation and Projection. J Open Source Softw \[Internet\]. 2018 Sept 2;3(29):861. Available from: http://joss.theoj.org/papers/10.21105/joss.00861](http://paperpile.com/b/ImpoB4/kwNI)

16\.	[McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction \[Internet\]. arXiv \[stat.ML\]. 2018\. Available from: http://dx.doi.org/10.48550/arXiv.1802.03426](http://paperpile.com/b/ImpoB4/udIE)

17\.	[van der Maaten L, Hinton G. Visualizing Data using t-SNE. Journal of Machine Learning Research \[Internet\]. 2008;9(86):2579–605. Available from: http://jmlr.org/papers/v9/vandermaaten08a.html](http://paperpile.com/b/ImpoB4/uvA7)

18\.	[Becht E, McInnes L, Healy J, Dutertre CA, Kwok IWH, Ng LG, et al. Dimensionality reduction for visualizing single-cell data using UMAP. Nat Biotechnol \[Internet\]. 2018 Dec 3; Available from: http://dx.doi.org/10.1038/nbt.4314](http://paperpile.com/b/ImpoB4/5APh)

19\.	[Ward JH Jr. Hierarchical grouping to optimize an objective function. J Am Stat Assoc \[Internet\]. 1963 Mar;58(301):236–44. Available from: http://www.tandfonline.com/doi/abs/10.1080/01621459.1963.10500845](http://paperpile.com/b/ImpoB4/RNVB)

20\.	[Apeltsin L, Morris JH, Babbitt PC, Ferrin TE. Improving the quality of protein similarity network clustering algorithms using the network edge weight distribution. Bioinformatics \[Internet\]. 2011 Feb 1;27(3):326–33. Available from: http://dx.doi.org/10.1093/bioinformatics/btq655](http://paperpile.com/b/ImpoB4/TzmI)

21\.	[Krishnan A, Zbilut JP, Tomita M, Giuliani A. Proteins as networks: usefulness of graph theory in protein science. Curr Protein Pept Sci \[Internet\]. 2008 Feb;9(1):28–38. Available from: http://dx.doi.org/10.2174/138920308783565705](http://paperpile.com/b/ImpoB4/xF0Y)

22\.	[Mistry D, Wise RP, Dickerson JA. DiffSLC: A graph centrality method to detect essential proteins of a protein-protein interaction network. PLoS One \[Internet\]. 2017 Nov 9;12(11):e0187091. Available from: http://dx.doi.org/10.1371/journal.pone.0187091](http://paperpile.com/b/ImpoB4/ZnMW)

23\.	[Xia C, Feng SH, Xia Y, Pan X, Shen HB. Fast protein structure comparison through effective representation learning with contrastive graph neural networks. PLoS Comput Biol \[Internet\]. 2022 Mar;18(3):e1009986. Available from: http://dx.doi.org/10.1371/journal.pcbi.1009986](http://paperpile.com/b/ImpoB4/9v33)

24\.	[Indriyanto I, Sumitra ID. Measuring the level of plagiarism of thesis using Vector Space Model and Cosine Similarity methods. IOP Conf Ser Mater Sci Eng \[Internet\]. 2019 Nov 1;662(2):022111. Available from: https://iopscience.iop.org/article/10.1088/1757-899X/662/2/022111](http://paperpile.com/b/ImpoB4/zIIm)

25\.	[Pal S, Chang M, Iriarte MF. Summary generation using natural language processing techniques and cosine similarity. In: Intelligent Systems Design and Applications \[Internet\]. Cham: Springer International Publishing; 2022\. p. 508–17. (Lecture Notes in Networks and Systems). Available from: https://link.springer.com/10.1007/978-3-030-96308-8\_47](http://paperpile.com/b/ImpoB4/gQAF)

26\.	[Nguyen HV, Bai L. Cosine similarity metric learning for face verification. In: Lecture Notes in Computer Science \[Internet\]. Berlin, Heidelberg: Springer Berlin Heidelberg; 2011\. p. 709–20. (Lecture Notes in Computer Science). Available from: http://link.springer.com/10.1007/978-3-642-19309-5\_55](http://paperpile.com/b/ImpoB4/hWVv)

27\.	[Steck H, Ekanadham C, Kallus N. Is cosine-similarity of embeddings really about similarity? In: Companion Proceedings of the ACM Web Conference 2024 \[Internet\]. New York, NY, USA: ACM; 2024\. p. 887–90. Available from: https://doi.org/10.1145/3589335.3651526](http://paperpile.com/b/ImpoB4/0Src)

28\.	[Guyon IM, Elisseeff A. An introduction to variable and feature selection. J Mach Learn Res \[Internet\]. 2003;3(null):1157–82. Available from: https://dl.acm.org/doi/10.5555/944919.944968](http://paperpile.com/b/ImpoB4/AGP1)

29\.	[Lee Rodgers J, Nicewander WA. Thirteen ways to look at the correlation coefficient. Am Stat \[Internet\]. 1988 Feb;42(1):59–66. Available from: http://www.tandfonline.com/doi/abs/10.1080/00031305.1988.10475524](http://paperpile.com/b/ImpoB4/xDZc)

30\.	[Silverman BW. Density estimation for statistics and data analysis \[Internet\]. Routledge; 2018\. Available from: https://www.taylorfrancis.com/books/9781351456173](http://paperpile.com/b/ImpoB4/KiI9)

31\.	[Hintze JL, Nelson RD. Violin plots: A box plot-density trace synergism. Am Stat \[Internet\]. 1998 May;52(2):181–4. Available from: http://www.tandfonline.com/doi/abs/10.1080/00031305.1998.10480559](http://paperpile.com/b/ImpoB4/GNeu)

32\.	[Chen YC. A tutorial on kernel density estimation and recent advances \[Internet\]. arXiv \[stat.ME\]. 2017\. Available from: http://dx.doi.org/10.48550/arXiv.1704.03924](http://paperpile.com/b/ImpoB4/K5Hn)

33\.	[Cohen J. Statistical power analysis for the behavioral sciences \[Internet\]. 2nd ed. London, England: Routledge; 2013\. 567 p. Available from: https://www.taylorfrancis.com/books/9781134742707](http://paperpile.com/b/ImpoB4/8xMS)

34\.	[Hedges LV. Distribution theory for Glass’s estimator of effect size and related estimators. J Educ Stat \[Internet\]. 1981 June;6(2):107–28. Available from: https://journals.sagepub.com/doi/10.3102/10769986006002107](http://paperpile.com/b/ImpoB4/80FL)

35\.	[Durlak JA. How to select, calculate, and interpret effect sizes. J Pediatr Psychol \[Internet\]. 2009 Oct;34(9):917–28. Available from: http://dx.doi.org/10.1093/jpepsy/jsp004](http://paperpile.com/b/ImpoB4/o3aH)

36\.	[Spearman C. The proof and measurement of association between two things. Am J Psychol \[Internet\]. 1904 Jan;15(1):72. Available from: https://www.jstor.org/stable/1412159?origin=crossref](http://paperpile.com/b/ImpoB4/g1sV)

37\.	[James G, Witten D, Hastie T, Tibshirani R. An introduction to statistical learning: With applications in R \[Internet\]. 2nd ed. New York, NY: Springer; 2021\. 607 p. (Springer Texts in Statistics). Available from: https://link.springer.com/10.1007/978-1-0716-1418-1](http://paperpile.com/b/ImpoB4/KacF)

[image1]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFgAAAAQCAYAAACSnrg1AAABqElEQVR4Xu2WvUsDQRDFt7bwo7SzsLBRLKy0thdE0E6wsLISrLWUiIggNkLAQkRE8APB6v4BO+1EEAT/AjsLndHd+G4yu3eJG0/D/uCxb97u3d5NLrkYk2iFdxkk4qI1mDOpRBtsyADYMamxPybUwPTkRuBVBgA3d1mGEl40YUemTlogjTZWdBc3pDHSPemMNGz8T+GtDAS+4xqcgr8z3wfwOA1zVZOZ5heKTyGeweNa33G+nNkz4flP6uB58QrU3U5hc4gDGQBlPtAcLS3+54yb4vstM689kHMycIROyL9XzAzpCCcsvNE+1Bd2nCKdW79O6rH+wY7MNWkIao3M5H8GQioDrzuGegm8o+hcvvlc7opNMaF5zPCCtLXotSzkOwXvoV0T0yvqAdKsyJBt03wOhrM+GTCr4F/siGC2C/6JNAg1bnoJ3uFrqnaxseE9+BuU2dF9s7S9tczhPiif2sLXDHlCV/NfIckkaQRqvlHHI/i/gLyvjuNrahnvwOwQ/JUd3yCrkkVSvwx/k5qo583XCwzRXoQn4NfAM1uirhLt4UhEJEqDPwCNb4zoK40HFQAAAABJRU5ErkJggg==>

[image2]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAH8AAAAPCAYAAADNnxDSAAACu0lEQVR4Xu2YT6hNURTG18iISGRkINLLgIGE2U0iJc/EiMGVMjBBKQMGJyMD/1JKeMmfMlJIStRLJuYmhpIiA8pEmbA/Z63nO9/Z574zePd43r2/Wu31fWufve8+976z93lm85NfKe6kuKGFEQf3BXFCCwuFKyle2d+FziVf1PgPmet7Mq8Y5uKGOXZXLIQ1NDLMxQ1z7C7YaENawyIrB4599oXrrribYtrKOdFeqpb/sMEGbwe52h4rx0O88zaINSuFaO6DHPdGibl3aYHQuVSDt5b3AfwLGU/vFa+Ric+4SgsBLvzpOTp+pppStIy27LTBC18rmmEdi2Sg14t3ylvuu130S291/OCwaBxUF5MO+t7qOKdFc55bA/PNW/VV46yjY2cF8jekuwRzv1fT8jfiKuVag57IeMxJb/tWX/9t0mDS/RzqqwbPvL1m1fpKyuHvI/01xSfSIDc2YP+caDCbnqGx0AGYe5uaVvocS6l23T1GNch5AP5W0Qq8i2pa/XPFD6oJ9FlHegflOq9qrFk9AA9vR6wPieY4QLUauQmaKFpGW3Jzz3bIiUWpx2BP/C5eoH1Vg5wH4L9WcwA6Tug2P2DoB+KBXD/Vek6oUHi7wqoXP6K8C/SDBzk/HpGoXSb/pnsgfL5+inLANd7vP5Cfmx/A76tpgw9czBZvz1q1xgdRbTd5G+iY2h/tcs9r9KzakQcrKO8CXUigPr7gAIerH6Sxv+IvfTV5fP1eygHXdP0Aj0n1GK2pZnQuhvVH0su85e+I0TFDn/E23pCCNVYeSmc4ZtU97TnlXfE4xUM1id0pnlj+VWqJldcH+PJ7pMFTK/vluJ9is+eFVffFHuVN3EtxXM0Gblnzvov1BUes/sc3LTrAY/2g5/tTnKdagPtzVM1/Dd7Dgf6ix4wATY+zMSNA1/9JHOP8BuI084vqdtPzAAAAAElFTkSuQmCC>

[image3]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAIQAAAAPCAYAAAAlM2e8AAAC8klEQVR4Xu2XTahMYRjHn0I2SslHymIWUixssLJRig2iWAg1IUmsZKPUbCQl2cjOwhV2IlGSWz42NlI2kpKSIqTkIwvep3mee/7n/z7njJnm3rlzm1+9nff/e95z3vOemTPnjEiZv0HbXRoxohf8Wk4Xuj6frgYPEb2sq+uLV0E/jtELVfNW+Yxz0sXgIWKZDHZdg5h7lcTzLpDYh/TrjphuvEttnOUUcV8Gc011zjMspe0fsKxCBx9jOQXovP6BRV9KzY+tf9uy3gGIvu/wvp6xMd+lPL9T5X5Zv2U5Oiaidb0rmRMS769Z16g8tfy7KGfgMdZAxoZ49tolqGXwzp1oSfui/U+rgk+Q/TfIdWPxC4K1KDtV/oltsX4V3FvwzdS+QmaiOdjhunzL61sPWZln3on2iVB/nnJIS2qKk8hd2/Lcnj9TZt5LXtPsd5lnhh3ndeS8z+NepvaCHMLjNb+GvDi17da/Z9ton9WB4+zHWWiZWSq55zyBFqqKS1j0mYbkc3fKTuTVXbf+w9QeQU1pSrHeV6ntLVXbaO0gS8nn44yckrzu835K7RbVlEUS74OcDhxmvUmuQXZ0zFjgQrRwkmXiDgugJfmjoarVoXPri5+z3JzzgTISeXRR/YfEHonqhyX3nBGtbQxcHVr3xxM6zpGL+gh7zXrDhPBgp8r3E55D8wXK+syM4H0vk4v6TfLOTdvukKLO++vd6Vw05zXGnT8K0CH68uhwXbP/2vnNqb9q49ZXbkh+np36UZ5An2lYnJ3aT3KTic7j/yL2WUY4I1qba/1NlhHPz0q2PG5+al8gX0ltp+T/ZPjYz1PbSg7x8fjO8DG1bZD5mFX5UMkW/qz18bEXrbkB3vHMfqA0pDihVqFL+MtSHftZAHtYGPrBrGRpbGGROMIisYsFcZSFcYCFsYFF4jgLIvpAo/OP2Mxi0Ohi3rAcUcscytEXYmiZUYuZAvSRg9dsRl2/tanNSm0FF0bUoi+zfyT+W9wz/wCkOyKrTNzAAgAAAABJRU5ErkJggg==>

[image4]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAOAAAAAQCAYAAAAPv3P4AAAE1UlEQVR4Xu2aW+htUxTGx5NLQnmRByVFyIOSHPJwXhBJeMEDnZBCpIRCnfnGg1JE8eC8SXIJ5Ul27vdw3AoRjsilFKIQ6ztrjbO//a0x12XvNf+nP/tXoz3GN+Z1rbXXXnPuZbZmTZuTVdhk/KPCXmTpsaTKPrC6Afh3VjZr4hv2lJqGvkHeZf1lxpJs+jbHkKxc/9vJT1b3A3uSdOd+m+e3Ntr5e7JlSDbvc2q0zacaDXa95LY2OuyexdSkvKnCUP629oQAtONUXIGoD6bEFxCUaHMMJfrPtQn9JBUrHlWhIdfOVCQr08eRKtj8SxaR06fkQhWGgsG9oKLVeu7ElQD9naDiBGzEwe+iRP87VWiILsJzJWYeUGFikrXHsyq59s6wOBdppfhLhSFEAzzRYr0kJfp7urLXVdxgcvMac7O5jPxHyFcusXZ/F0isfKnChCRrjwccUNnhKmY4VeKoPQe59yXeSJbqTyudFmjgW6svaLCNdIfraH2NAbREflRf62l8ldVrWKBtuBbxGfl3k18CHcNHlZ1Dseb7joOWV5B/sPFf5kSGvvZWIVm7fZ0LL3MObjTwuc2vD6zxnB3kK3wNaL8bwfcq9PGz1QNNjcF/nvLON5X9SvE+5AOebKJYPx1d710sscPa25W9S/H+1q7TFwPW4L9FsTIbaF1wf4dJDG4STf1XKHatC+TdhjC03DIka89H6cprfGNlp4imoM6zKk5M7vg+VNl+KnaBRp4ItGsCDZzV+Gme2o0PKPcMrIONYn0kwxeD7+BRnfMoxpfzB4qB1gE+1k81UQgeg/etuHYo+a4/RrFrXVxr/WWYrrI+XpzXyHLzcZK156NA851LzWs8kzgCdfBLWhodG0g232UeRNQItC8CrY/LLT4h31X2jGhaRmOg2tgYf6nglzviFovHWgLuA/4uih0tE/lOpDFj5zWm7FiSDZvPw42PNezHjY9lAmLm8coOFE2J+ihBdG3dYfX+yWCiwUK7NdByIHc1xT+SD7zusYGmcdfJ8tgnnsu/JjEDbV+Ju0gDrQvuw/9zZY4KtC76yiJ/iIod9LW3Csm6zynoyzNHV3azisR9lm8Dm1fI3d7EqfnECwlYEr1U2SeNBlAWTyQz0cA20ph3VOjid4sHC+0g8vnT4Z04zeXi3IGG7zFvmefKsMa+x75m6us3ikugfQyJk9XrnbMXU7vR8kpfXhlbfgz32mL7uEDxUoCDZYeenx1Wz19//Zyu8SL3p4oNp1u7r+PJB1dKzP5Pgabk9BbYrXmjsg8r+01yX1vdkDbmWrRr6LmvNGFxWwCab4DgcUPL4G7HdaN2WHuRfEdj4HVmmigA1qS4s+pxwdo2mg/wscFwnhDrRpG+7QFwQ8Xj26tWP/bjrt7HFhv3azkGvOCB6wtr8z9I940o2BWkg1+svhZ8/tEx0hjgB2GX1dcA1va5LyHXZT/6q0MfJXN1mZy+ZpOQO4Gqa7wsU7UzBcnqfQQFf0EcQXEifyw83+cCzYleK+Ny71V2KcUObi5rNjHRxQBUP0biZdFfoL2NPikAnTuItD7OtMX/jJ2oLdeuC7SLrH4K9L0GJ2pnzSbkNqtPpltayM5Z9YSvWr8UPHddHq1Cyfnm1qtr/uNEa8Eh8Js4/xdKfgE7+Rfmnc8D8OB7ogAAAABJRU5ErkJggg==>

[image5]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJEAAAATCAYAAACQoO/wAAACqklEQVR4Xu2Zv4sVMRDHp7Pz7ASxOLDwB1gI/gGpzkILCzubhxZqJ1japBSs5Cy0slIQGxUs5QoRrCxEEKuzsbAQBBERBN2RhDf73Uwy79i8jbAfCJt8Jj9m53KPxx3RTG3+oJiZ2QufUczMrMr8abQGTqAQ3EXRGD9QJGjhEq1Sx82unQV3BMbNcQZFoIXiRxak56P5CP+SlOaMxXHSz9I8os1jfwVlC2gJ30YxEZzfTnhquTJHUQTip2xu7Rhwjm+75il/lqWu2vprpMcm5TWKQM1kHQoDpUuUivEFuhn6D7t2WsRq4SmdSyQXYy6jAErr184HFIKXKEbEoTCw6iU63LXf4HBODTzlzynVNbeW4fgnlExcKDc4BuMaaPs/QyG4FJ689o3o7wt9Cw6FgdIlWnTtFMoJ8JTPk8nVt7SWGczx4emoH+T+YzGW7NKyqKWWQ4tbPPdvJbwFh8KA5X1eoJgAT+U8tfgjFAra+n+BizCujXaG5iWWOUz84ZdaCcu8UjwHr31XaDxnIy5Q8FTOQ4trXpKtGUoc10A7Q/MSyxwNh8KAWjjBTxQT4KmcpxZ/ggKQ65J7SHkexshX6t/IXMuhxb+jAO6R4YUyOBQGSu9zsms3UE6Ap3yeTKq+X1AAuCeP+Z0HUvbviHEtMLHIBRp+SX1Fy/n83BGx66JvwaEwULpEuS+r68RTPk+uK9cXya3h2H5w34LvweJ+eA6CleA/pW+hDKRy+Ej9i/Sga++XYTMORQa+rLEmsiEpt04cDXNM5ZRyzAJFgOc/RRno7YX/u9IOqoF2lubHwKEYgV8oGiVV15RbGbkJ9w+IcW22UQjOoWiU5ygaRavnaJeIL84om+0B7VztXyKtcRVFo6TqebBrh1D+r2gXaRdFY3gUjaLVUau7mb9UeOTenTBhvAAAAABJRU5ErkJggg==>

[image6]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAOsAAAATCAYAAABxy/qhAAAD4UlEQVR4Xu2bT6gPURTHTzYiRD0sKRayIllQ6rcgGwvJn5IQxZZQlDTZoERiI+ntKKRkQUlvIX/CUlkgK1kiUZS45/3uePf3nXPund+8O68703zq9JvzPWfuOec3M+/Nb97vEXW0lc0oDMFfFDo6OuqhZ19XULUL75mx2Sh2dHTEZSkNXqBVLlam6n4dHf/hk6iHokU7wTS9yaxCQaHq7FX3i8kNFDxsp7TvBmai4HAFBUuP0jgOldGa1/TdKDQQbbb3KACvjS1HsSS3aXKfe4dhBukzajqi5R1Hgfq5uT2FWF3k9TS0mKY3Aq35ZShYtPzU6Rl7a7d9M7xDwbIOfP7Nc8RjLg/tq69uDMaMPTeWkb8W9ich7V9GC11EMQnVuYkChfdJGql5SWP2oNBQtPkYX4wZ9taQHy7lhNaORUb+Wr4YcxUFyx/wN1FxrV2CVhdc5zSKDlIfktYYpOYljdH0puGbQzr4nO9aWd6Av8HYNdDqICN/nz9RAKR9rxubBtoc6ueOONpKq00Fl8hf6zEKJOTnAgbQTwGpp1soWKTcquyzr7zmC2d7ut2uk9Ac31BoGBn5Z1xs7DCKDh9QIP96LpxXNjcGoVr3wFfz3YDvJ05GE0OGLDa45jFjB0DLwdyquOvw9llBr5NQnVA8dTIKz/ADBYv2kC20Xk7ZvFhwvScoOmA/6I/D9/14Uh51/NispX4NvvXyGTaL/hgN96ccvj0K1eX413wHQFpTgvPKmPZ5yyVUMxSvi9D7mL+X/LTXR0bhGbT4sLpLKGcRlZ+xDJzHD/p8+RhDfxwU0U8F7OuBsfWg5WBuDOpYM0SoZiieOhmFZ9DiF1CwaPk5oXhs3Hq+2hhDfxwU0Xc5Q/14GYsNrslfDshvSxHMnSwHaXDNV852nYTmCMVTJ6PwDPdRIP8+vthv8HeCHxvuJQNf6+87+GKeK+ItcUpIfWkPWKTcYeHPF/k6/DrmxPY723Uxi8JzXEahYZwj/4z8cGkERfLvo8XuUP8Y5vaZ9NwY8Nq4/kJBY/g50VbQpDy6aOwRTSx+fjCcDFLzksbgT6mqfKLBC3aU+p9V6iY/Fq4hv1BoEBkV55NmlLR5xragCJwCP6NiLa1mDJaQvjbrqwUNKWhzwS8kJITUm+8hzQkUWob0frQNaUZJQ8rkpITUb0FjIX9I89LYXSeWGoXmLcPqbWA+Ci1E+m4vU+a47kUhYTaiYCnMycIaY1+MHYJYahSat/hO3FEUWoL2veA2cRIFww5jC1BU0M6X1OB/nJBoSv8i3Lw2wDYUHPgfsduE9h60iVjfTBs2f6r5iIJl4Fz/B3+GW39rCpCfAAAAAElFTkSuQmCC>

[image7]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGEAAAAQCAYAAAAGeBHHAAABrElEQVR4XtWYvy9EQRDHpxKNWnsFCaIj0V6lUlH5ByT8CaK4hkQvEZUOjULjRyLRSSh1Og3+BCVv7nbO+L7ZtS9uz9wn2bydz87tztu3uV9EPnhD4ZQlFAk+UXhmH4VjxlEkeEfhmUsUzsk54YcoPHOLYgSYR2GQ86DcMFLFKn6r+wxFLqvUm1wWWFBjpYjdDPvlcGUeqrZCf7i5TLAejIWYZ15RKGR/sfW5UH3+sNxQMYKTxFon5Md4QlGxo/oH9F0kX0/V2KA5J9gQIxZinomNsR8L/fUQJ9lGUYA22Q/pXvW50BMVl4TXmlPxTHAWMc9soiA733J9plEUok32Q9AkCx0wuBbH1+AEzBWaeMtRC2J9KhCeIKd1Qn4M6+1IYxZaCFwLY01szPJbVPdHhuueSm6aWlIBrDXEram+9iW4ofr8En/8sD0wl1ms2izKAOZj3OU4XHlQ2jCw1hF3pfqP4VoKvGf+FibxrvJCqm4LGZtUfTfwCfRAk42J/VhrMoc7PPxt0WQDrdw9FKPGf/+Bd0f2xlpMoAjkvt41zyiGSAtFgikUgRcUTfgCR0SBiSaII4wAAAAASUVORK5CYII=>

[image8]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEgAAAAJCAYAAACR6WqOAAABFElEQVR4Xu2UMW7CQBBFp6OCjiIFkosUNBwgdQ6RggKJMqdwC5QIcYZwidSpaFBSQBspKdJQI8GOPIu/PxMDkt35SaPZNzPy2ivbIj7HK141G8n3+C9vLTPan9p6hQ3Jek8hZrZGEqillnnGpSuXg+zIULL+LeExsMwz3+R9cqXMdZ2Q81rzo1MvRYdewBdWqxvew/MJ+XuIjxAHyQ408mp9hF3hGrsLD7HXBe9zryPa+wWfh/gCVz5DrMF1puyaZ3goOterZBfiDbwtxf164HvL3v3Ef5D2WlCPs2Onxv5QqBL8auo6+gjqVePd7BL8WbJPPYWaoj/3CF/jxzI+Az8bEv2vUCV0qMPFhhw+1QaiOSDiBPYtZPlxn42CAAAAAElFTkSuQmCC>

[image9]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAKQAAAALCAYAAADx/yZ6AAADjUlEQVR4Xu2YS6hNURjH10QmkkQy0fUYSHnllVC3EFIyEZLyKANvKUmpM2CAGJGRRJEUIkTRCRmQiYgBGYi8Jp6hFN//7vWd+9//s87Lta9y/eprrf//+/Za+3X2XvuE8B/lZ4xjmvhD+PgrNZHAa4dqoguU1egiF9ToCeDi4cJ0F0XP1cr4rdQ2w++Ol9puRUj7/zw46GlqFkjRJ7mV8VupLYpnFufU7Ml050XpFYqdb31ofvw9ofnaIvmr+4DJU1E0mGOkxTKLEnkcz6PvwJtscd7inkU/yj2MeYB2Vmxxwynwp1icjX19Guy3uGQxx6It5M8H+mtCts/oD5Y8gB5h8TT2N+XTOZDvH1vE2Hw6HLf4EPs8z/vQuZ2D/hjqb7OY2pmugNyo2PL2n8hDfKbcouilgD/a4qTFdvJLMQfQ4pzXGqMDPZjlpFPwzjaKenylvh7oKtGOeqwPkMf+INGgkT4T8j8E5H+QBvOiD9qpD3Q81YzmUpo/tjZT3+FtXlrMtNgXNX5UCCY1h6Kea/WBelqLlmsOW7wlXUEH6p3wisJ3crgmQubvEK8UfYZ1OeGB+yF7cjq7QnVNK7pvbPHESn2VXw31t2fwhNdcPY019bvYx48YPLB4EvuM34Q63uOEp/pmyN4cTJ/Yai300oSHp2+ZNANdEi+5VjmU8IoEc3mor8C7Qnp89BT1Uho3jYMbjGuavaFSHtCTjVd+vVo8aZ0t0XPuik6h+ROiNQ+NfXKwNEnVpLgRqnOqgXqNdAcwdye8CeIpqGk2ajFXNGoninb4VcxAHxQPr6pUHSiJdr5YXCet+85vDd8XoOM46kPvFM9J1c4Q/Z10Ch2D13x4mh4lDbQeekDCc3h73RaoN068W6KB69y6fa/FBjZC9YZFgXmw46ydYaId9VyvE69Eemv0ABbroNY4WGwDnMBHsQ+Q95q22OL1peM46rvemHMzatXiwwAsDvkbDPir2uEx8IHIeG5SwlPtc7KnuM95rU3pdtIfQ3aOgT7NKxufCul1SFHgteg3ANZ3WG8w2C98XPHrzH18uaK9bDEwZOswzitYC/EaEqAOT+TXIRtntsUbyXP7jfoA2x0hzQyxeGUx3eK2xYuQHd81LoqgFhdnfsi+pDEHvtjRd/yC1luiLLBYHdLnEd5p8nCs5ZCfE+DV7cBbaHGRPPfxI2ewfrwTsn87au0fs8RiLfu/AIx2YWHssZY7AAAAAElFTkSuQmCC>

[image10]: <data:image/png;base64,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>

[image11]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAQcAAAAQCAYAAAAbD7owAAAFBUlEQVR4Xu2bW8guUxjH1w2SY4gcU0QJEcWFi69sCRcOcYGtiNhxQ0kOFyYiJedDCe0SOeSctkPqUxIXyiFJOewLhwsiXEqx/mY9vv/893reWfN+M/Pu7/3eXz3Ns/5rrVlrPWtm3pk184awoI3NKixo8I8KLXylQuIPFRaMypIKygehnmw2w9PnGW+cGouPnLxl0vtmv7BtPwz0J6cPxScqODyuQsLro47vQifvAtLXGjpGjoWnD8GxKnjkOnNFtONFm2c+ViEDYvSmaL9KemhycwVy2lBcqYJDrk8ldVHv5Iw2T+Tm8aRo54s2JNp+lmtDs+AR0fag9HqgJFA6oc+TPybowzeSHpuSNh9TIZTV0ziX1Flr4EKg4zpX0kODu4e9VcyBjr6X/I2krwcOiXaGig6IE+6obtSMEfkprBxYeoCNRVu7X6uQWFbBwfb/ZUOdL/gi+CRnjEjbPP7Hn6Eu+IxmrAOKApTQX7VZMet+3BRtTxWJXN8+VWECsx7fGLwb6jH+oBkjUhxjFDxLxZ7BI0sXDlaBsAOoxCbRlq90LT8EJeMaEtyOPqwicZEKoVt/q9Ct/FoFY7xbxREpinHpibRauu6/a/lp6NIGJnKMOE0Cbe+ctn1wVLTXVSzAa7+rruB5HI9uKL9F8uaJvs+5afbTWscKYBGytfAq6br/ruWnobSNt9L21VBep2+43T4PrGnw2u6qK2en7azHNyQ6j7NiYtua2TYhG6K9omLkiWiniPaapME9kq6i3Sva/tF+CfWHGic2sxpYX0tsEm35QMsg/ahoQ6N9eCSjGXtFOy35t0c7IfmI9eXJN/g1bhXtyOR73ykA77Hi8GgHqJjw+srcJmnUOUY04/60xeu/O5N/XrQHk89cHW2TipGX0xbfTpzKGQOicZh0jL6QtruG5l3U2+SDO0J9B2gspy3m41nSFa/dbAZO8pwO9iXfyuwe7T7RLgkrX2Dxvt4gH3De1rR9MdRvD4DXj75payeXP2lCh8BrC/r7KkZ2Cs068O2isEz6QSmP0Xo5vAVJrzxoW5D8QoXgx/n7UI/j95ReCs1yfAFsG49puby+ybWBt4M5/fS09fqvdSxdpa3FBvxGPqP7+D/gML7aHCd5WhGdhcavXbQM8AbAPu4MqlA/63rv7HP7HgL80i2J9mPw44A+c17uoO6LSe1wnvYReLHUspPSmmd01Q39iExjyegcVI1cv58PkL9PtJ8p/R35xl8q9EwV2sfhxeDA0PwidCv5XPawaNdRGh8x7pJ8xOAuymO0vVXzdLRzkp/buTdpns+UlBmCMdsaC+9AYh+//rdQGvDt93PkM7l46eNKjly9afHGxP7n5PPJY1wTtu8P/ngs/Lh0WagfGQyNqxcPBh9B5e7+OoOO2NoA1gQMbrgSDc9LV1EaB+u3yX8nbQ37VbTPkfHxy6XBH1jflHw+vZbAOoOxY6hjCTB3uI2tUtrie2va/p22APMFPiQNeJ9Al8yVV3cacGdheCcDfkEN07GOZX5Jn2eJNy7zscbCaSNXVvH0qcAfYTarGOp1gg2U3iHazcnHGobxEPngzGg3iHY0aU9xxgj0GqwZw78yS+QDfT1Ykf8S+cAW+hjvj1efqeAwdJzxuhlrYcb1YeXxmRdZsaALcKzuRvr2BJ8zeu7xnVoV7VBJGxeH5vIBKP7j1YIVdAIWNPFO7K4X8sVftmfLNm+A/gX+AKPE4LHAXQAAAABJRU5ErkJggg==>

[image12]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAAALCAYAAADP9otxAAABTUlEQVR4XsWWu0oDQRSGz2NYBMFGUucJLMQyprCwD6QLki5os62tIDY+gAEvhSkT8gJCEFL5CDYhCMFKnYMz+OffMzsGstkPDsx858zOJTtLRPJ8+3h3ce0i8xF82Ry4eJK/+TLITcH3vLsFd+mdcg++Bf4GvImVKBxQEtacd9QPaF2bpeTHB2Jeaiwcb1IwoGRw3gtoM1o3ILfv/R75M+qvkFG/L9VtXglvwSEnCK2Zkxt63yGv1/pf1KXazSsn8ruGBicIvi5j8I/gv6CdZJ3N48cpFeuiY2YsCX72MfgleOuKm+jAU5YVcOXiVdIHhweAtTFfiBY+kzui/jZ4gbauqegajMTeYDiAHU7E4FcpYDlk01eA61Jju2LXfHr3Qd7EeoCi7oFliVhrUGJe2RU7H/4kJVnI6q/EsQ1wPjzwjHKx9UxYOJouzlkiP+5JiggwLQhmAAAAAElFTkSuQmCC>

[image13]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAPAAAAAPCAYAAADakUJeAAAFkklEQVR4Xu2aWch2UxTHV1IiIhEu9H34Qi7MEhd6yixToQylT0oypRBx4RSRIT4yFr0hQ8oYZbh4L5QhU0Qi+Uy5MJTMijj/96zV83/+Z53h1Xve9+nz/Gq31/rvdfbeZ59z9jl7P4/ZJPeW6esyzZfpszL9TWV/lOmDMr1fpr9I/z/wj6ehWG/Dt6FEe33bXEzsjBUGF+pjFW26LuBy92Xo9lD/U4m2h2hLCeq/WMWEOPchxiC7zzYUflChB0syUf5q9UrUX0m2suXvz9DtDV1/Rt82+8YtlhdsuLqngRU9NzR+KNnTxIdlelfFAXnChh+DoevP6NPmcdYv7r+AevdUcQNhfxtu3HqBxiOtFNH+Nu6vs2p9Dg05kvK8NfcZ+iluX+o+1v0ZD1tVfoXnm08WL9A0PufbWN/H7bPGxQuc6TrSXp4Hx4vPHGCTZS+RDY6x5n5tb+OyMzxvY96qrzHkSBmvW1XPTVpQcrpVZYXoz1rzdcQei/Zrb/FR/rjbt7uvtI1DgLJj3Y77gc+jz/FafqFV5/Oz53xusCM+xn+3cfGCfzP5wcs22Q7sgvxGtHNdXGJVxX1SF7iQgfZD/YD1uLGCsJFnOsPaGvFBIRrs+8jf1qob4xf3n7R6fCE+r393JJ052qrzigv6zGRxLZ7Xl+dZvQ99178bqVjyoE3WtwXZAGV4iNhXVFvvOfR9Sec4TAqhhc52+Ey2zsZED7JjY5JnjRmJBvtT8kPLeNPGkw+3/SVpjLYTE4zG1YjKOwMHADvgBfnvkQ2yPkF7lfxdXAtu8VyPVX8u0RbjX0Xai6QHI2s//nfP/yzTd6R3gTreIv9Aq97+QVubbTTFqf452ZdbvVz9yxItxkt19sPWmEDHYZVNjgNzjdXrgY+vruAI15guH6j2k+fQ460fbOL5TlY/jmkrmyAC7yR7uUG7kZgHyvSJaEDj4Osbpu/FwFs3eMy1IOuT+iDTAPSnE03JtDZ4vL5oKFOtD1kct4XEXx8AGj4Rg7WuMfD10zjg2IPEDzINcL90HJSIU03918jPPvHVB5kGmnSAsmtVdHpv3OITjcFB2MTpw2lW/1RuSn04zOqDnJ3EDlbX1QfQCvL3c43J/N3Fj7c5WO2akmlA9a0TbbNECw5WwWmKByi7h/yjXOsi1swKtGytFugx8DF5qpaByZnvQcRhncuc4HpGk56B2CLRuvx3Eo152+o/CQYay7SV4b8X+C9GK1erYFWlbRUvNVl7GJCAy8LeWHQQfhbP/m1u30oaE/795GM2DOC/4XboTTc+UB1+XGx89od2t9sneQ74raZoveAGz7UMvn6dZHxbpu9VtOp4TDzKnOdZe+CcRAP8goDOG4YRt7NoV5LPaNsgxkHR2DsSLfxvyD/RbYB78yO3Q+c6fiQbk6jWH4xsXMYxq0nDznbwCNkLICCrfGSVrgv7odA+NPm6q8txfC4Xic6E/2iiha2Dihs/3sCxzsaDe6NrANpD5DNZ/dicKkQH3K/tbLxxxCnQc8O/5gL0OTZ/Drdxn7Eb3AbisCGn4FcBba+pL9xP1VUD2PXGlxEoLI/TYxgt43Fg1lo9Fv4q8vG185Xb2KAFh5TpFbcB3sb41YCfD653U7Kh8+Ycg0n8XNEQH5PDb6JPNRgQ7LTyTMfUZh8Hb9E1buMfTOuoDJwtPphXoeQCq3ZtA/0c2tUmH1Csy/mtXJCdcXKZricfb3d9ULRffIPgQo/ID3B+GBueqQP0eY589LmLrhsFu6nXqejgbVuQr5MFfv+NDT/lVBt/gQAdC70eCsaqaRwCfqMHPD7BXWU6UjQs2WInGeCnSL1+z4kPsvoZ/HSpoJ0t3S5scvk2Y0bKyPOuB3jGjBlTBn67xoOLJcGQ/8OeMWPGQOABxj/FZkw5/wKR4QPLdogjQgAAAABJRU5ErkJggg==>

[image14]: <data:image/png;base64,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>

[image15]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAGEAAAAQCAYAAAAGeBHHAAABrElEQVR4XtWYvy9EQRDHpxKNWnsFCaIj0V6lUlH5ByT8CaK4hkQvEZUOjULjRyLRSSh1Og3+BCVv7nbO+L7ZtS9uz9wn2bydz87tztu3uV9EPnhD4ZQlFAk+UXhmH4VjxlEkeEfhmUsUzsk54YcoPHOLYgSYR2GQ86DcMFLFKn6r+wxFLqvUm1wWWFBjpYjdDPvlcGUeqrZCf7i5TLAejIWYZ15RKGR/sfW5UH3+sNxQMYKTxFon5Md4QlGxo/oH9F0kX0/V2KA5J9gQIxZinomNsR8L/fUQJ9lGUYA22Q/pXvW50BMVl4TXmlPxTHAWMc9soiA733J9plEUok32Q9AkCx0wuBbH1+AEzBWaeMtRC2J9KhCeIKd1Qn4M6+1IYxZaCFwLY01szPJbVPdHhuueSm6aWlIBrDXEram+9iW4ofr8En/8sD0wl1ms2izKAOZj3OU4XHlQ2jCw1hF3pfqP4VoKvGf+FibxrvJCqm4LGZtUfTfwCfRAk42J/VhrMoc7PPxt0WQDrdw9FKPGf/+Bd0f2xlpMoAjkvt41zyiGSAtFgikUgRcUTfgCR0SBiSaII4wAAAAASUVORK5CYII=>

[image16]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAALQAAAAdCAYAAAAD3kYUAAAD10lEQVR4Xu2bTahNURTHd2JAQhlIGZhIyVcZSBm8DAwUAyEDI4kykDIy8W4kYuRrIGQmDCgDBgYvAwozZYaBiSIkAzFhr3fOdpf/3eucfc7e5+M+61ert9d/r/151nvn7HPvM0ZRFKWEhdZuWfut1ropDaAbq8woNKGVGcV2FBRlXPmJgqKMM/xxAw8ssaYorbMR/BQJmaIPRamMlHCpEjK2vSJwF4WW6PsFlea32gwT+jLUVUUagxMSM+6cYeWg9Up/UXxam3Q9fhFHUWB8N/KepqTp/ttmPQo5fJ1UXs58Edyc0+B3wTIUesIHFDy4hMZ9TclSFGYAvv06zMpUv4T5ItgR+l3Rl3lwQufUZFI30Wcf2IYCELTub9aegzYJflcELaBl7qAgcMgME3or1MXSx31JhbQ2SZ/GVfr+gqDPoS/iEBRDB8a5ebkq1IZee2Fb9I9b2wlalzxGoQS3v7iuGPZa24RiDr+uO6ytMd0c7AdmdM3oS/jiuLaHlafBBmW+4zorU8xsVq7CKVbmbY+A78C7BzIVaCGcQwHwza8MauMsBT9QyOH9U/msR2+Dr/lPHBd9CYzj+4d108JD5q/NNQ76PkJiyrhhRi+C7+1BirFCWGWysegthcRFFAIRL0gNQvooiymqd/P8VGBU/8U1EOBjrAOfQN8h6V4wmPyTHq0Ier1SFhMC9oG+oyjBUjGwtsLafSPP4zUKFVhpsn5nYUUNpPlxQmKahB5NPzOf5nOB+UVUmjsGO3+RR5Og+nfg1wHboU9ssHYJRYDahVgoUqykh0Btj6FYE/qlK8IdRh0vWZmYAL8JcK/Q323keWBsITx4M/O5/oqVHTwpiga8aorrOTyODi6+dlMm+6+QNqF5HEDRjCZGKL+svUcxgsUmu6Nw6HmZXx/aNwd/jLtmssMiPV41CV5L7rvyTaZxsG0pPDkpkfjiHahNmNGkJtvvAnImcj2EgcneGri+TvxTmxHaV0r4Oh23wQ/lihntKwW+PvH6vLX2bFj9F1/bJqBxnprsrQyOiT7nCQopKBowBfzTH0Ia7yMKLVD3oIz4+klFTL8xbUOZx8o0Ht1VONIc7qGQkoMoBCJNlsNjpHhJbwMaewHzH7FyCPtM/Py3oADE9L8LhYTQvF7k5Ulrb1idY2DtAYom/gtdhfieI1MxZbKFSxeFPrDpEj63OrdAaV1VCOljDgo9we1fle+bhKw3GjqN/o+cN8MNrrrRVeOR+Sa+j3FjOQpKelxSlX2AwIlJRPeuuujOpSi1oaSq8skgT8ZYU5TkVEkuTMhYU5RG0OTqKX8ApKbK3sOPdG4AAAAASUVORK5CYII=>

[image17]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAARQAAAAQCAYAAADXloDbAAAFXUlEQVR4Xu2bXahtUxiGR7lCSiTJTwflL0TcyM0oSVKSUiIikgsX3OBGZyWp0ym5kUS5EEJRyF9qXaD8XcmNUtt/RCiFCx3Ga42v/a53ft9cY+2959r77DOfGs3vfcc3xhxjrrnmmmvMtVIaGdk6JmqsmBfVWMC/alR+EM15UZutYhrESg7irSYHsadXSt/BWRVramyS30TjZBvihMulXK3miMuPavTwtxqFfWqk8YLixeB80SvjzDT8C9HCVo/B68/zNgv6nIh3YfV3CjtpLJ+q4fC4GpU31EjjBcWLjaGPxyHNsWl1Bxj7uVvNbeIqNbaZltfAy3lXjcp4QenGxrlpdt6PDABOtufVHIihT+xl2EljAWeVcpmawgtqpHge4wWlGzPhMUHFnaVcXsrXTh2+x99Wyu9SB1B/Yin3ij9J3R1C49YS25NLeaWU46v2gH9LKU+UcuR81RzIw1cBjJMX6ODfV8ol5IEbax2w8Zt+K80O4JOlfFw9AzmvlnKO43vAP65uPyjlcKr7vPoA26NL+aZq2y/3i1iL8ZxooHoIon18mWbn07Vp/jyI8ltonaPnGd+rUYna6PFv4e26tfx7Stlb4z6mQazkII6YlHJl6o5ftZKDmHH7YHO/aG0AfbNo42XSumWuS11fNVAPJ6kH3qjGUWm93bN1i4vW+zUGPLaPyP+klJ9JAx7DtG6R8xD5QMcK1GO9lzweT18bcJrjcXtGNWP7ail9LKoHloMPq5Z8D2un7VUDzzOius/UqHB+1Jb5lmLk/0XxIqZBrOQgjpjUrY5BtZKDmOn0AaNjVrw69RB7C2En1C3n4kJintevass3HeH1B3AHALTu1rpVH1rXJjjH7ny0HVAP+iXSuXqGxdoO2LxxXLX+pzT7pGb21K3mqh4Cbx9ror0cj748e6qgOaqB5xn8YWjk1F34NrzXrJVl86dBrOQg7uOYND+e10v5g7RHDmKmM8eOQXh18HSy8KwwerdjwHvU8fr0IqIxAM/DVxb1Pc13Nobm/VLKe+JpzjJjM7w6zwM3pPk6fG37k/RQROMx+K51EYvyWucY9RP5Oe3+CwrGwksS0IeR9shBzHTm2DEIr469mygGmq/aUB/aFsr2kNcCPrl4bQXt+MWAfq3Gdtdk/tOkzevTAN4DjqeoB/2UeBdVP8LqdE3IeIxi+FgXYn0EaQX1raWPlvqzSV9B8bK0zjEaU+SDrfrKY1yQ5vN1zc1jGsRKDuI+dOyqPXIQM51+OkZa97QO+gvRTKT5K9Gl5BuqwTLehPQ7FANrc+qc2+0LC7fqmeYfTJl3huNx+6iva8R7hDSD7+JrNX6GfO3XUF/1UHj76TsOm0H7Um14/sWlnK4m4bUBy8wF9ZaD7a9St4hpECs5iPvQ/av2yEHMdPrB6u93Ncab/R+qw1MTezNhO1mv+h/uTH85CqyeX0h4D5I2j7cA6wi22Hp7mi2semCdxJ6cYP2D33zA+rT1FEMPBLR+elrOKY7n3TXcQZ4t7uY0vwjN+9UxMPen2d2X5pj+cM6dzd2OEXK03VB4++G53lXK9Wl2J6EX+2VpmSPONX2iB7xcJqpvfb0A6jHPaZrdNfNxaAHtvFjJQdwHxoCnqSfV+MB8tUsOYqZ1biMHMXiR31RzIPQivCqiOUYneOQb9qhXWeaCslmGvKAwrfPIQWzgh234qcPILgMnyCpPfGUV+2udo/fE8WE1AnbrT+8xbrujxjpi6zxyEBut/YwcZPSta+0WWubo/egSRPnKbv1z4EbnkIMYbNufA0cODXj9aDvgJ17KV2r0wBcusNE340YY6oKyUXIQg/NEp/8A9Of9X7MDlVUAAAAASUVORK5CYII=>

[image18]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAOUAAAAVCAYAAAC5WykPAAAEtUlEQVR4Xu2bPagdRRTHR0EEC9EUEsGPp6QQK0kkoCg8MIqiiEUKIYWQSklSRMukuIUWYiQBwY9ClBQB8aMzwUIeBIWEgEoKbQxBkBQhTQgRDYG4/9w57577zzkzO3d2H/de9gfD7vmfM2dmdmfeftx9IQzMOzc6LAMDAx1Qu6ieCPU5BgYGiK4WVW39gTnnGgsZRixsILWTsbZ+F8ii/IsdhfzIQoeMWFBsYmEJKZ0npfE3K9wVt3eS71ey2/INCxtAauDHWHD4goWWnAjp9sFBFhzuDZOF+Qn5umQtTNrB/qryvaB8h5QOvid7mXiABQc+16OovUe6huskkWCuxHYJn7PQM5+xECaTSiZdjvMsFHBbyB+vP1hI8FuY9L1PkN8a96gpK6QJe1iI9N3XjaDNGL5iIdKmrjVPTbxknt6W2volpNpquyhTOXKg7s8sGnzHQgJZlDX9ymHltjTB8z3NwgLzHAuEdww8XdMmJpwMduBpFogtan+r2tdYefvg5aa8waKizaJ8pin3sRimn49wW7lZ2RoZ6/1NeV07iNJj0ufCfClM590d8s+DXj88HdxN9uNk98EK2W1vS0FqLOAtFsL4tlXq7dUOAvMU8zUJEuFkMF7H9GJFzH9qn6l9UdEWq20N/GssEl6OS015tSkfKM2KhaZ1KwZ4egrJ/Sc7KtF95v57nGUh4tUVHQvx97h/Sul9cD5u8QglfwDajg+k4jwf9ANk71C2xsuxjhfg6XeofR3zk9oXcFDwwiDFqGVJ4fVVgH+NRYJzHI1b62Sy/aKhsS1gkZeCk231oxbkwzglL7bbJu5bwLn0rgJW346ofXl5BQ6H/q6Wuh9fKzt3/jXWWATPB13fZcDGOC28HOt4AZ6uycWMYumbXD/gT52U91lQcO6PDA02vyXlGCHVjxTIJ6UrkOsVslP5R2H6Da0mVQ/gLX4uBuhx5koOHh9zkYWIlxu3nY+wGOE6bGtSvpt4AZ4uXG/KvywSi3KlTNVnH9uANdjeFWWWK6XA7dSgr5AaaA+zGCm9Umpy/j6YtU2vnqefacrfpHmxIOULXwY/wNNX4hZ+3JII1u+Ai/JM+S4LkU/DrbnFvmBobN8zpY7h2Lag3rMsVoB8Vl88XSh5ptyl9i1/H+xT+7pN/bwneH3y9F9YiCBen+tzTflB+RhLWwfOVRYjayyEcbz81vYP6Rae3jWrTdnJoiI10TwdwPemsh+KGngnbvHcwDnEZh1YWg7U+ZbFSpBTv7wS8ByY6qPns3TRMEm134rtCn3sdTsjtS+6tVCB1b/UHQ7Hi/38lDoG83SVRQ0nY3L+HLX1S5i1rVnrzYr3w7MH+pf6SmSj8Y7XkyzMOd44gHVHkoovwc0jDjcgkvOnmIcvenK8zULPlHzRAy435SqLc8AyfNHj9dXS8Tv8dhZnxJynx8O4Yf1MlGJZvn21KI2vpe23r+DDUN+/2voey/DtK47NYywG+yODro5jMs8KCxl2sJAh2XjPfMxCgiss9EhJW4+G+mOI3wH3s9ghXv/kWXsR8S4k3lhLKJmXA3PGa6F+EqB+bY6BgYFI28V0e1MebMpTYfJvY1wGBgYq4UVVUwYWgP8Bgneh4iKCKhIAAAAASUVORK5CYII=>

[image19]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJAAAAAQCAYAAAD59vZgAAACs0lEQVR4Xu2Yz8sNYRTHz5KUhY0FykrITlJKvUX5kY0Ff4ESpSzsWMzCxko2shA7SlK2yLvyhthZWJCVyA4bpfAcc073e7/veebHe++MeSefOs0533PmeZ4z89w7c6/IsPjNwoho01ub2sGzKP00dImFEdK0x9ssrGYOSD8b6BELI6RJj31c617Rhi6yOGfuszBi6nod5Qbqmj7mGApVva5NdoTF1cZSssJ8bbaq4XlRNQfnOO4TnfsgxCdMa0NVfS53L9k+809jYmiclekmFihmCplssir7YvU5frJg7JHl83PcJzy3xhdIq+MTCwCP76Ceq1kpfo+acFVqajW5geK2F2glLLJg6PwFxJtNQ7gmh18o3cyRNbmQnOf4JPicc3K96qNrDYuGr+0UJ+bEMxYqyPX1F05y3BWvWDB4fo3165wpWOgAfXRE68nFb5Mdg9h5zoLBYzG+ierq2nKOhRoq5+ekx6w726T8RNVZdNOR3Pise3x4Su1nA+nc1yHG958boDua28KiLO/JaapjfCXZVylfNfQb7CbkvifbKNP1+6W8Vi9B4/HfSzneN9A+J1ufbJeUT4EsOJj6Hl8DvQu4CWWTTOs7Ieb6guIu0DnPU3wHYobX6ES6boQceLOV3CZWoutzFHx9h3Gi2sh/HGiVaNFH8/XY6KQZieZQ7ZYdPV+AjxQsdICvQ+0BaUykOVEu0pzLMpkHN7CD57qP2odAw/gXaXshdi3yB4c+DpE2iy1YmDP6LtN0PT/A5/90uEfnCQstwHU9teM70KJN9dqOx5Otk0l91CNvoO0QDwpePMc5FqSsfWN+F+j4+pdCHVqHxkTaCxZasCPZQ/P5RrN2xo76YdBvNWW3HR2vP5RsK2l6zl2IBwm+BA5pofNYC/aGzDL2LOf+ZyREP/Wb8s830B9TBuZkHrsWrwAAAABJRU5ErkJggg==>

[image20]: <data:image/png;base64,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>

[image21]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAANAAAAAQCAYAAABnTPHAAAADd0lEQVR4Xu2aQegNURTGb0lZELJANqKIWFGyUC8LIQsbslEslJIiJbGRlLKysEFJssCGsrJAYWOBslUWSkSxkNhxv/+7970z35w7d+7735n3Ju9Xp5nznXPvuXfuvJl5854xZXosdIC/LDhCepcIzeEjC1Mmg3ksOM6y0CAvTP/EecQBhWMsOEInnicWbxuMJzSmbyw4QvlTxoS2IH5htVgT3Bf7y028rhbfxIKC1m4cPLH2R/jauPawILjGwpTxsMTaBhYF2sLmpmfKdeBvJU3ylgVT7oPBFT2W0xY8jkOKBjQNhPT/mbEck1jRWDwXW8ivqvuYBYe8ojM/3Laq37ZYYfRx1NXAJRamBI9Vo8SKxuJNgJovWRRoY7rCAjHXbbntG2u7rL23tt3aHWurTTkvN1r/mnbO2iIWHVp+CORijkutfbJ2xtoNpzXNNmsrzXC8OOYgZfx1SOkPucdNf73xqIxHYjzxpPQxQ6xBLJ4T1IJ95wChfViqxiljnHfQbVlnPzfc/wFFA7gzX2XRoeVXge+Z3IZ9xq9JHQvhY9ieJH2V8GdL1RiYNW7LbdiPEmsQizdFqO5sdfbBYVPW2c/NeVOscZp8SarOXHZb5ONu65njtLbgWuyDp0bX65Da7oO136RxH5vdNnh3UkVBLN4UqKvV1jSg6adYMHoetL3Cn++0Nrhthm/bQjVT9RCcD7/qUTkn94xeXyOkx0hth/xlwj/qNIn0sb9A+AOxilg8B6jBdTQNV0x5oks4N4SWxxr8m6Q1zQlTHgfI+QjH+exr+HWoY1UgLt/2Xnca2Cd0EOsrRGo7zoeP72uSz2If8fXCH4hVxOI5QI1fivZQ0UJo34s0tD5Y8/7CgpoX1HhAvkbKSwT2JRdNOe79LwW1GbTa69w+fveTcG5dUttxvvd3FNQhnD+DKlqemX7M26tiODtyHL4mo2kSfp6V8HzgA1z9uF/vs54Tnu9O4UtCY8BbNCZ03AD03Yq21pR/QmgCHlfVMda0OqS0w/r/JM23f11Q+wT7xtWt6ofUSWG/tcUsEsFJdpjQnFL1LjHqHEZtF+Or2D8i9gc0VTgndcZY5688XWKUv/JcYKGD1Fnrtrhl+m9MYXdlgCm9XZgw6h7UunldIPXPpM9Z6CA9snHTIwuykYUJovKTrxA6wbpEaA7vWJjSPv8AeVEfOLXiZwYAAAAASUVORK5CYII=>

[image22]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAMwAAAAdCAYAAADmQIN4AAAECklEQVR4Xu2bzasPURjHz85GioWXf0BC2NgokZBEiayUrpKklHeSNAs7pLAgFnZEkaSIiAUWwsJCUhYoL0kkZeFlvn7nuI/vPOfM/PjN/OZez6dOM8/nOWfOmXPn/dd1zjAMw+gJP6wUimGojM7LWpaGYejY1dQwusBOGMPoguMsDMPQ+cLCMIw49jhmGF2wl4UxdLGrX72UzW/GooWU7YNkG4u289B1dvAJJyLsZNEnXrjOuDPyMda4+I9fA3lZL+KDTq/XBKl+R7DoI6n5BDdYRIi1byUZxWWDL8s3xSGx/tqVj4vziFeKOPMulJsi1zSXWAh4P/oFj4PnE2QUx7jDoq3MdfqOs5Okck2BMexQ3CxyEh73R3KZWO8nT1kIlrLoIzyfbxUHNCe5wmKogR1cwNIzLy+TWPYB7aRGfJmcBPkJIt7nXSAT6/2E90uSyjUNz+di7xjNScryUdBwVV4OuM7jRpaXmd43xRmX7k/LXXWDHstpfjn+d41m0MaWAvW1Eyb4hYOpJN32WwYeL2NofV33S+Qm5+WsiJuE5zNwjQWxheL3eRmTl3d52ZSXk65zYmrbdstdMcFxXcx2nb7Oc0LwkuJ7fol2H4THwdbUuAH6esSyBB7fEopjB0AAOXz8kHOQqh9I1ZFzqKG13eCXnOO4bmL9rcvLLpYerU24a3GO419AnlZcjNuuk8dExwry30ODCqD+RZYu/QWNx4j4GTkJ8lXG/SY0SHDEFfsvo0r91a5aPdR5xTICLkapbaZyIJbf7Io5jiXhwpz6G4R8FVL15rr4B5RYOxx/X8mpdVniBfAzubqJTZTmApzjuC5m5OUcyxJiY2M/UnGS7a6TP+GK70MaIZ+qd58FEWsLL79QTfWuCcr6WZSXCyxd8YuaBNscJeKj3hVgyXGvwfa5D82Bwyw8u12xfoi1ieolPOnzKWZ4nHfFOuc2Kk4S7iqfhHss1mNgm7dYuvhVWBIbD3vEe/w63gPqgvuV8xnAO/kUlq7YVsI5xGP9Oi6Sv5jjE5IQ4xZZB9wfgMOLPLsYyMkDVT52aBPYK9AHFwlivDwGnnvHJaC1Z6dRpY4ktl3NMbE67EM8IGWPKZvPgOYA2sfgNiE+JCU2gF+XJbFB9BI8BuDKiEcL9IU7BpMag5are9z8R9L602KtSBDv90vMR11wv0BzjFZnoit6vD+i8B24l/A8avMJNHeMhWDADX75C2AbuGmMI99K8G+yK1ga/wQOgN+PFq7z2FKVZSxajvZYrp1Ew4ZhvXN9gq/G3cyxfGdqO7H9wo+cw5ZvLIx/Zrr782Cq+lk6EDsQ20TGwjMUxm60kHDgnJLyP2ArC8OoQupF2TAM4oGzE8YwusJOFsPoAjth/oKf0oqpIRshzsAAAAAASUVORK5CYII=>

[image23]: <data:image/png;base64,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>

[image24]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAXAAAAAQCAYAAAAPmiTRAAAGvUlEQVR4Xu2ca6g1UxjHVxQlInxR3jouKVKuhUg7cr/kklKi7d4r5S5CbRIhPrgUUSe3REkRKbSlpBdJlHvOFyIld1HE/N9Zz3mf859nrVlr9p45w5lfrWae//Os296znzN7rdnHuYGBgT6zJQuJ/MOCh/UfyO6KP1loiZ9Z6JAN6vwMdZ4Cv09JTF1Z0SrPq7g6FtymevPi/KJMWJwjn7rqnHUJ8RgLBOqO/DHEj66c21ekDww0TeBHseDsZB27Ltukq3676sfiabIvIjvGDWRPXDkX5KlarEnXJTJmvcuLF0J1oF/K4pw5x8X7Z6bO1gX49vXHh8gnSH0cc97ggbVBkwQeuiZ3ZqHgfRY64jMWWuJDFjqEEzgIvTeaP8jmOrD3I22ZRVetAL5xth4CsWMWE8jpY96g71D/MV+IlPiUmIG1y7wSuKWB/1ICD80hRt8S+AMsGPA82f7J0JYJJaqQHiInVtO03jxA35ez6MmdP0iJT4kZWLvkJvCTirIFiy58nQ0JvD2sBA5iY9qsKIeRxvO+w9CWgeNE0j7wOgPt96Ks8+fsY65y5Z38Ia66LiwJUgrWvMEzyq+B/VJRtiVNc4HX9jJ8FrEYGRdAe7srm+G5WHHsl3kKMnZ87X1K6TsUZSfvA5e46jedk13p53m/V5TjivK5K8f/RFF2oxgN9KOLcp6rbsLABy1nb2TeTFx17Gx3Dfo/19CakpvArb62K8r1LHpSEzi3i/f+O9Jy6FsCP8ZV22c7l1ACj7Ub8wmICcaJUxcL6DsqW2+QHOmq9Z4lzWob9og0Qce+qY5TpeuYm8gG25OtwQXJ8Ro93o+UFgKbmzE/gP9RFp099nf8UXQcdyUNSOLXyLwv9Ef2sw20pucutkbGZiF1U0ouqHOZoa0m3P8bRfmVtBxyE/hfLBTc78o//BZNEzjsEWk59C2Bo+0rDG0W2kzgJg+6qnPR0M4k7WxXTaQ3Kls0tk8wNM2CP77lVvre9UeO13bMZwF/LAa+Uw0tBHycXJhQfdbZBpYGWGcbd4esWfatpAnwjcheLbjvV4ryN2n3uWpcm3BfsDcnbdEfOdYiJ4E/zIIn1k9KAj/eVdtg+9iifGLoIfqYwDV46IA1sXHcXzsC5CbwI1gwQN2tWBTg5MatO1PYX5Km4XisK7PGNrA0kKqLfZo6TwXxofXv7121PdjXkabheAsrJmXs2AixYlLqwi934lqL2ZqYr2t4LGwLUxZa4kpXHYNl70F2jJwEHmorpIOUBI76jxuaRUhn+p7AYWOpV/OIOud4i9wEHtKFOv/GgAMMjSvCxjprCImXdd3flCaIPfHHi5XGiH6X0rCezvEH+mMowcWIxVs+0bCGbGHVYawYfFAsXQP/QSy6tHmzH/ZzysZryDGamM8C8aklh9tctQ7bwpSFlkD/2GBiDYz9EXsWWJMWQmMWchL4ayx4Zl1C4TE+qbRttMNVY0P0KYFjX4jbFntrpd2uzjneIjeB/8KC4kWy7yR7I1bD0FiHzQnkHnUu8W/7I3y6Ddy9i73oj7DxzDTQP2hZ8gXou4AvXPkhFr5W54DHDO5lwWN9yxDq9Dp/DCsGf0AtXWuWX7B8+hl09rMNYlrM1yXoU/e7oGwez5TstkC/OuEueQ3wtx4AHydAJjWBf8uCYtZNTH49YV9LmsCxIfqUwPlaEk0fNdjUT7mmchL4zSwocHPLrGgDa8oyCTwloi8qbBRak5Fz/GXnAcHGj2JYA3iiRey7/TlAkjnUVdeNsUwRSrCi6TtIAb6D/TnmwOvXgszbKh+rOAb+0JLL1c4er+ZlVz5FY4G6e/pz/LHjccTa5nmPN7nM11FsrF8K2Hgb+3PsZej+saYs36zY1yXyHoFxUXZRNjNloSX0mJaUxk8YgX18qSM1gYfmLoT8dQkcTyGhrtx94hyf1bE/ZyzNoq8JfOTKX7GKvbc/alLHEUrg1qZ/rE0ZH5eBlsCLezqLxPAGzEbO6zdloQXwfueMSdjAApGawPHkUozQ2OoSOOrpb9d1hPph+pbAU4ktczBWAg/9kOcFFga6Ry6E2AWREjMQ5xaX9/p1kcAxHqxvpzDxJWXPIiWB17UhNPkpfWrbQmp8kwTehLoEjic/UseMb1ITX15d4bGxErjV1+ssDKwO8ubw7rUGMdYHaSAdvIZnsRhgRKUtrA+mxcgoMeaZwJv8M6s6vzCiUkdqu7NS1w/8Kf+DaGSUOjiBh55aqxvjwMD/ilNY6AHXsDAn6hI49koOZzECJwsrqWv4F7jzoi//TrbNa0kvj2EvIQT2zBrxL8bQeCPRXotyAAAAAElFTkSuQmCC>