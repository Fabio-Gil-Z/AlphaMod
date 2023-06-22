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


## AlphaMod Pipeline
![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Modalpha_Pipeline.png)<br /> <br />

## Results_
![AlphaMod Pipeline](https://github.com/Fabio-Gil-Z/AlphaMod/blob/main/Images/Table_1.png)<br /> <br />

