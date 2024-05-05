# GenBio_Final_Project
**By Sarah Couture**

## Study Background 
The relationship between reproductive and systemic physiological function has long been of interest to biologists, though the ways in which function is coordinated remains elusive. One potential mechanism of coordination is the co-regulation of genes important in different functional domains. One specific example, the potential co-regulation of gonadal steroid hormones and vasopressin (AVP) may link reproduction with the maintenance of water and solute balance, both critically important processes in animals. In rats, neuronal expression of vasopressin has been shown to be estrogen dependent, though the relationship between non-neuronal AVP and estrogen is currently unknown. In addition, it has been shown that an increase in non-neuronal AVP gene expression, important for water retention when access to water is limited, decreases sperm counts and motility in male lab mice, and negatively effects embryo development and litter size in female lab mice. To highlight the potential connections between dehydration and fertility in mice accustomed to limited access to drinking water, we have conducted an RNAseq-style gene expression study using a model desert-adapted mammal, the cactus mouse (Peromyscus eremicus), which can survive without water. Using specific nuclei important to the maintenance of water balance (SON and PVN of the hypothalamus) as well as the pituitary, we quantify differences in gene expression in groups of mice exposed to two different water treatments – one where water is freely available and one where water is withheld completely. We use these data to as a first step towards understanding the impacts of dehydration on fertility in a desert-adapted animal. 


## Methods 
<br>-Extracted (n=20) brains from female (n=10) or male (n=10) Peromyscus eremicus that were either hydrated (free access to water) or dehydrated (no access to water for 72 hours). 
<br>-Extracted the RNA from each tissue sample (SON, PVN, & pituitary) utilizing a typical Trizol protocol. 
<br>-Performed library preparation to isolate mRNA of tissue samples and amplify each tissue sample. 
<br>-Performed bulk RNA sequencing on all tissue samples to observe the gene expression profile of each transcript. 
<br>-Indexed and mapped samples utilizing the UNH Premise cluster.  

[![](./RNAseq.png)](https://github.com/SarahCouture3599/GenBio_Final_Project/blob/main/RNAseq.png)
Figure 1. The bulk RNA sequencing pipeline 

Raw reads then went through quality control and multiqc platforms, followed by the selection of particular samples following a PCA plot. Samples were then mapped to a reference genome utilizing the STAR index, and all visual analysis was compelted using the R platform. The code can be found here: https://github.com/SarahCouture3599/GenBio_Final_Project/blob/main/Code.sh



## Results 



[![](./sampletable.png)](https://github.com/SarahCouture3599/GenBio_Final_Project/blob/main/sampletable.png)
Figure 2. Sample Table showing the head (10) of each individual sample name, the associated gene.count file, the specific tissue extracted, the treatment received (yes water or no water), and sex. 

[![](./PCA.png)]([https://github.com/SarahCouture3599/GenBio_Final_Project/blob/main/PCA.png]
Figure 3. Figure 4. Principal component analysis of gene expression of the pituitary (pit) and hypothalamus (sonpvn) of P. eremicus. The axes are labeled with the proportion of data explained by principal comonent 1 (tissue type) and 2 (water treatment).  

[![](./plot%231.png)]([https://github.com/SarahCouture3599/GenBio_Final_Project/blob/main/plot%231.png]
Figure 4. MA plot showing the differnetial gene expression of pituitary and hypothalamus (sonpvn) in both males and females with or without access to water. The log fold change being on the y axis, and the mean normalized counts on the x axis. 

[![](./MA%20pit%20df.png)]([https://github.com/SarahCouture3599/GenBio_Final_Project/blob/main/MA%20pit%20df.png]
Figure 5. MA plot showing the differential gene expression of the pituitary in both males and females with or without access to water. The log fold change being on the y axis, and the mean normalized counts on the x axis. 

[![](./MA%20SONPVN%20df.png)]([https://github.com/SarahCouture3599/GenBio_Final_Project/blob/main/MA%20SONPVN%20df.png]
Figure 6. MA plot showing the differential gene expression of the hypothalamus (sonpvn) in both males and females with or without access to water. The log fold change being on the y axis, and the mean normalized counts on the x axis. 

[![](./Bean%20Plot%20Pit.png)]([https://github.com/SarahCouture3599/GenBio_Final_Project/blob/main/Bean%20Plot%20Pit.png]
Figure 7. Bean plot showing the log TPM on y axis, and sex and according treatment (F - Females, M - males, hydrated - with water, dehydrated - no water) on the x axis of MUP4 gene expression in the pituitary. 

[![](./Bean%20Plot%20SONPVN.png)]([https://github.com/SarahCouture3599/GenBio_Final_Project/blob/main/Bean%20Plot%20SONPVN.png]
Figure 8. Bean plot showing the log TPM on y axis, and sex and according treatment (F - Females, M - males, hydrated - with water, dehydrated - no water) on the x axis of MUP4 gene expression in the hypothalamus (sonpvn).

[![](./pheatmap.png)]([https://github.com/SarahCouture3599/GenBio_Final_Project/blob/main/pheatmap.png]
Figure 9. Pheatpmap of the differential gene expression in the pituitary and hypothalamus (sonpvn) for both hydrated (yes) and dehydrated (no) samples. Expression level (0-20) is on the y axis, and female samples are the on the x axis. 


## Discussion 


## Future Directions 


## Acknowledgments 


## References 

