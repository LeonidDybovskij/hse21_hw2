[This README in russian](README.ru.md)

The purpose of this work was annotation of the genome of a bacterium using the related strain T. oleivorans MIL-1 (http://www.ncbi.nlm.nih.gov/nuccore/HF680312).


To do this, GeneMarkS-2 and E-utilities were installed.  
Also the required for T. oleivorans MIL-1 data were downloaded:

![Screenshot from 2021-12-19 13-53-57](https://user-images.githubusercontent.com/60808642/146672975-a2b3917b-cacc-4442-98bb-6e26363a9c46.png)

BLAST was also installed: 

![Screenshot from 2021-12-19 13-54-10](https://user-images.githubusercontent.com/60808642/146672984-c506ca25-111f-49e3-98ae-24f35820aa8d.png)

The SwissProt protein database was downloaded too.  
Seqkit package is also required:

![Screenshot from 2021-12-19 13-54-22](https://user-images.githubusercontent.com/60808642/146672989-fba50718-9713-4a82-9492-b3a13b24534f.png)

Then the work began with the scaffolds.fasta file from HW 1:

![Screenshot from 2021-12-19 14-07-25](https://user-images.githubusercontent.com/60808642/146673001-b907ddae-d51f-42f3-b45e-5d307be60d47.png)

After receiving all the necessary files, the processing was performed in the .py file (src folder).  
Unfortunately, creation of a variable for each feature was not successful, so due to confusion around one variable in the loop, the necessary features were not added.
