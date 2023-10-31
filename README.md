<p align="center">
    <img src="https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/DISAMIS_UNISA_tri_logo.png" width=60% height=60%>
</p>

<h1 align="center">
AlphaMod<br>
An Automated Pipeline Integrating AlphaFold 2 and Modeller for Protein Structure Prediction 
</h1>


# Abstract
<p align="justify">The ability to predict a protein's three-dimensional conformation represents a crucial starting point for investigating evolutionary connections with other members of the corresponding protein family, examining interactions with other proteins, and potentially utilizing this knowledge for the purpose of rational drug design. In this work, we evaluated the feasibility of improving AlphaFold2’s three-dimensional protein predictions by developing a novel pipeline (AlphaMod) that incorporates AlphaFold2 with MODELLER, a template-based modeling program. Additionally, our tool can drive a comprehensive quality assessment of the tertiary protein structure by incorporating and comparing a set of different quality assessment tools. The   outcomes of selected tools are combined into a composite score (BORDASCORE) that exhibits a meaningful correlation with GDT_TS and facilitates the selection of optimal models in the absence of a reference structure. To validate AlphaMod's results, we conducted evaluations using two distinct datasets summing up to 72 targets, previously used to independently assess AlphaFold2's performance. The generated models underwent evaluation through two methods: i) averaging the GDT_TS scores across all produced structures for a single target sequence, and ii) a pairwise comparison of the best structures generated by AlphaFold2 and AlphaMod. The latter, within the unsupervised setups, shows a rising accuracy of approximately 34% over AlphaFold2. While, when considering the supervised setup, AlphaMod surpasses AlphaFold2 in 18% of the instances. Finally, there is an 11% correspondence in outcomes between the diverse methodologies. Consequently, AlphaMod’s best-predicted tertiary structures in several cases exhibited a significant improvement in the accuracy of the predictions with respect to the best models obtained by AlphaFold2. This pipeline paves the way for the integration of additional data and AI-based algorithms to further improve the reliability of the predictions.</p>

# AlphaMod - [Paper](https://google.com)
Cite this paper <br>
<pre><code>-->put the citation reference here<--</code></pre>



<h1 align="center">
Figure 1 - AlphaMod Pipeline
</h1>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Figure_01.png)<br /> <br />
# HIGHLIGHTS
![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/highlights.png)<br /> <br />

<h1 align="center">
RESULTS
</h1>
<h2 align="center">
Table 1: Average RMSD for the different types of models obtained (domains only) Test set B. 
</h2>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Table_1.png)<br /> <br />

<h2 align="center">
Table 2: Pairwise comparison of the top-ranked predicted targets across different methodologies including: AF2, OP1, OP2 and OP3.  
</h2>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Table_2.png)<br /> <br />


<h1 align="center">
Figure 2 - TARGET T1038-D1 - TEST SET A
</h1>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Figure_02.PNG)<br><br>



<h1 align="center">
Figure 2 
</h1>

<h2 align="center">
CORRELATION BETWEEN GDT_TS AND QUALITY ASSESSMENT MEASURES ON BEST PREDICTED STRUCTURES
</h2>


![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Figure_03.png)<br> <br>






<h2 align="center">

<span>A:</span> GDT_TS vs <a href="https://alphafold.ebi.ac.uk/faq"> pLDDT</a> &nbsp;&nbsp; <span>B:</span> GDT_TS vs <a href="https://swissmodel.expasy.org/qmean/">  QMEAN</a>



<span>C:</span> GDT_TS vs <a href="https://prosa.services.came.sbg.ac.at/prosa.php">Prosa(z-score)</a> &nbsp;&nbsp; <span>D:</span> GDT_TS vs <a href="https://saves.mbi.ucla.edu/">PROCHECK</a>


<span>E:</span> GDT_TS vs <a href="https://salilab.org/modeller/"> DOPESCORE</a> &nbsp;&nbsp;  <span>F:</span> GDT_TS vs <a href="http://molprobity.biochem.duke.edu/">MOLPROBITY</a> 
</h2>

<h1 align="center">
SUPPLEMENTARY MATERIAL
</h1>
<h2 align="center">
BEST METHOD DISTRIBUTION ACROSS DIFFERENT EVALUATION UNITS
</h2>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Supplementary_Table_Evaluation_Units.png)<br /> <br />

<h2 align="center">
Supplementary Table 1 - Average GDT_TS scores for the different types of models obtained (domains only) for Test set A.
</h2>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Supplementary_Table_1.png)<br /> <br />

<h2 align="center">
Supplementary Table 2 - Average GDT_TS scores for the different types of models obtained (entire proteins.) Test set A.
</h2>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Supplementary_Table_2.png)<br /> <br />

# CONCLUSIONS
<ul>
    <li><p align="justify">While AF2 has achieved remarkable accuracy in predicting protein structure, our study has highlighted the potential for further improvement. We have demonstrated that, in principle, by combining this cutting-edge deep learning tool with traditional modeling strategies, it is possible to achieve a substantial improvement in the quality of a protein’s tertiary structure, especially in terms of GDT_TS. Only where AF2 fails to achieve high quality results on average and top-two best comparisons over these targets: T1029-D1, T1043-D1, T1047s1-D1, 7MSW-D1, our AlphaMod procedures cannot significantly improve prediction accuracy.</p></li>
    <li><p align="justify">Furthermore, as described in Section 2 and Supplementary File  4, Tables 1-4, large-scale protein predictions can be effectively applied, thanks to the automation integrated into the AlphaMod pipeline, spanning from data retrieval to automatic processing. Finally, our pipeline provides a unified platform for comprehensive protein structural quality assessment, encompassing several metrics. This addresses the current challenge where these tools are dispersed across multiple service providers. AlphaMod, on the other hand, offers an integrated solution by centralizing all these quality assessment tools within a single, easily accessible platform.</p></li>
    <li><p align="justify">The current pipeline is only the first brick for the development of a tool that will also handle heterogeneous information, in addition to sequence-related features, to perform better predictions for selected subsets of proteins, with non-common structural features. According to our research, the addition of supplementary data has the potential to improve the predictive accuracy in most of the predicted models.</p></li>
    <li><p align="justify">Moreover, in future research it would be of great interest to study the feasibility of jointly using supplementary data and AI-based integration models to improve predictions in situations where AF2’s performance level is below 50%. </p></li>
</ul>

# FUNDING
<p align="justify">This work was supported by UNIVERSITY OF SALERNO, [grant numbers ORSA208455, ORSA219407, and ORSA229241]; by MIUR, [grant FFABR2017 and PRIN 2017 program, grant number: 2017483NH8]; and by BANCA D’ITALIA (NDA, AM).</p>

# ACKNOWLEDGMENT
<p align="justify">This work used MARCONI100 (https://www.hpc.cineca.it/hardware/marconi100) at HPC@CINECA, which is provided under ELIXIR, the research infrastructure for life-science data.</p>


<h1>CONTACT</h1>
<span align="justify">That would be it for now, if you have any question / suggestions feel free to send me an email to:</span> <br>
<ul>
    <li><a href="mailto:fhgil@utp.edu.co">fhgil@utp.edu.co</li>
    <li><a href="mailto:fgilzuluaga@unisa.it">fgilzuluaga@unisa.it</li>
</ul>
<span>*Thank you for reading, have a great day!*</span>
