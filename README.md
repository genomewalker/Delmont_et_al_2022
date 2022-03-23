# eggnog-smag

Code for the some of the analyses in **Functional repertoire convergence of distantly related eukaryotic plankton lineages abundant in the sunlit ocean** from Delmont et al.
To recreate the analyses follow the instructions below:

  > This has been tested on R 4.1.1

Clone the repo with:

  ```bash
git clone https://github.com/genomewalker/Delmont_et_al_2022.git
cd eggnog-smags
```

Then let's install the packages we used to plot the figures. First start R to get renv installed:

```
R
```

> If you open the project file `Delmont_et_al_2022.Rproj` in Rstudio it will perform the same steps.

If everything went well, [renv](https://rstudio.github.io/renv/articles/renv.html) will be installed and you will get a message like:

```
* Installing renv 0.15.4 ... Done!
Successfully installed and loaded renv 0.15.4.
* Project '~/Desktop/repos/Delmont_et_al_2022' loaded. [renv 0.15.4]
```

And restore the environment:

```r
renv::restore()
q()
```

You will need to download the data from [here](https://doi.org/10.6084/m9.figshare.19403531) and decompress:
```bash
curl -L https://figshare.com/ndownloader/files/34481117 -o data.tar.gz
tar -zxvf data.tar.gz
```

Then you can run the code in the folder `scripts`:

- `01-smags-eggnog-lower.R` will process the outputs from eggNOG-mapper so it can be integrated into the AGNOSTOS results
- `02-smags-eggnog-get-coms.R` will produce the tables used for the analyses in the manuscript and that can be imported into anvi'o
