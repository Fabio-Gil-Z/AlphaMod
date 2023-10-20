<p align="center">
    <img src="https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/DISAMIS_UNISA_tri_logo.png" width=60% height=60%>
</p>

<h1 align="center">
AlphaMod
</h1>


# Abstract
**Motivation:** The ability to predict the three-dimensional conformation of a protein serves as a key entry point for
investigating evolutionary connections with other members of the corresponding protein family, examining interactions
with other proteins, and potentially utilizing this knowledge for the purpose of rational drug design. Considering that,
during the 14th round of the Critical Assessment of Structure Prediction (CASP14), AlphaFold2’s 3D protein structure
predictions demonstrated a GDT TS score surpassing 90%, signifying it attained an accuracy level comparable to
experimental results. In addition, one year later, the team responsible for AlphaFold2 released predictions for over
200 million structures. Despite achieving remarkable success, only a relatively modest proportion of these structures
(approximately 35%) were deemed to be highly accurate. Therefore, our investigation was designed to evaluate the
feasibility of improving AlphaFold2’s predictions as a post-processing step by merging AlphaFold2 with the traditional
spatial restraints approach applied by MODELLER, resulting in the development of a novel pipeline we called AlphaMod.
Furthermore, a validation of AlphaMod using two distinct datasets was conducted. The first dataset consisted of CASP14
targets, while the second dataset adhered to the same criteria as a study conducted by another research group. Finally, this
study aims to introduce a novel bioinformatics pipeline that seeks to improve the accuracy of protein structure prediction,
with three primary objectives: (1) identifying protein families that are best suited for the AlphaMod approach, (2) confirm
the feasibility of improving AlphaFold2’s predictions and (3) offering a comprehensive metric assessment tool for large-
scale protein structure evaluation.<br><br>
**Results:** The coupling of two distinct modeling strategies in our pipeline has enabled us to demonstrate and validate
an enhancement in the GDT TS score in two different datasets. Additionally, when conducting pairwise comparisons,
AlphaMod yielded a superior GDT TS score compared to AlphaFold2 alone in 70.4% of the cases in the CASP14 dataset
and 56% of the cases in the second dataset

# AlphaMod - [Paper](https://google.com)
Cite this paper <br>
-->put the citation reference here<--



<h1 align="center">
Figure 1 - AlphaMod Pipeline
</h1>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Figure_01.png)<br /> <br />

<h1 align="center">
RESULTS
</h1>
<h2 align="center">
Table 1 - SINGLE DOMAINS - TEST SET A (AlphaFold2 vs AlphaMod)
</h2>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Table_1.png)<br /> <br />




<h1 align="center">
RESULTS
</h1>
<h2 align="center">
Table 2 - WHOLE DOMAINS - TEST SET A (AlphaFold2 vs AlphaMod)
</h2>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Table_2.png)<br /> <br />


<h1 align="center">
FIGURE 2
</h1>
<h2 align="center">
CORRELATION BETWEEN GDT_TS AND QUALITY ASSESSMENT MEASURES
</h2>


<h2 align="center">
<span>A:</span> GDT_TS vs <a href="http://molprobity.biochem.duke.edu/">MOLPROBITY</a> &nbsp;&nbsp; <span>B:</span> GDT_TS vs <a href="https://alphafold.ebi.ac.uk/faq"> pLDDT</a> 
</h2>

<h2 align="center">
<span>C:</span> GDT_TS vs <a href="https://salilab.org/modeller/"> DOPESCORE</a> &nbsp;&nbsp; <span>D:</span> GDT_TS vs <a href="https://swissmodel.expasy.org/qmean/">  QMEAN</a>
</h2>

<h2 align="center">
<span>E:</span> GDT_TS vs <a href="https://prosa.services.came.sbg.ac.at/prosa.php">Prosa(z-score)</a> &nbsp;&nbsp; <span>F:</span> GDT_TS vs <a href="https://saves.mbi.ucla.edu/">PROCHECK</a>
</h2>



![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Figure_02.PNG)<br /> <br />

<h1 align="center">
Figure 3 - TEST SET A - TARGET T1038-D1
</h1>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Figure_03.png)<br /> <br />


<h2 align="center">
FIGURE 4 - OP1 BEST MODELS CORRELATION
</h2>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Figure_04.PNG)<br /> <br />


<h1 align="center">
RESULTS
</h1>
<h2 align="center">
Supplementary Table 1 - Average GDT_TS scores for the different types of models obtained (domains only) for Test set A.
</h2>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Supplementary_Table_1.png)<br /> <br />


<h1 align="center">
RESULTS
</h1>
<h2 align="center">
Supplementary Table 2 - Average GDT_TS scores for the different types of models obtained (entire proteins.) Test set A.
</h2>

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Supplementary_Table_2.png)<br /> <br />

# CONCLUSIONS
<p>Despite AF2’s outstanding achievements in terms of the
accuracy of protein structure predictions, we have demonstrated
that, in principle, coupling this innovative deep learning-
based tool with traditional modelling strategies can improve
the results obtained. Furthermore, our pipeline provides the
capability to benchmark a wide range of existing datasets,
due to its robust capacity for handling large-scale information.
Unfortunately, where AF2 fails to achieve high quality results
(GDT TS <50), the use of MODELLER cannot significantly
improve prediction accuracy.</p>
<p>However, the current pipeline is only the first brick for
the development of a tool that will also handle heterogeneous
information, in addition to sequence-related features, to
perform better predictions for selected subsets of proteins, with
non-common structural features. According to our research, the
addition of supplementary data has the potential to improve the
predictive accuracy in most of the forecasted models.</p>
<p>Moreover, in future research it would be of great interest
to study the feasibility of jointly using supplementary data
and AI-based integration models to improve predictions in
situations where AF2’s performance level is below 50%.</p>

# HIGHLIGHTS

![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/highlights.png)<br /> <br />



*That would be it for now, if you have any question / suggestions feel free to send me an email to: fhgil@utp.edu.co, fgilzuluaga@unisa.it* <br />
*Thank you for reading, have a great day!*
