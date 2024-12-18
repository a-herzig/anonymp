# Anonymp

install the required R packages

``` R
install.packages(c("data.table", "dqrng", "tidyverse", "matrixStats"))
```

download 1-User's dataset

```
mkdir anonymp/1-user/res
cd anonymp/1-user/res
curl -O https://lysine.univ-brest.fr/anonymp/data/chip.txt
curl -O https://lysine.univ-brest.fr/anonymp/data/coordinates.I51.2_15.txt
curl -O https://lysine.univ-brest.fr/anonymp/data/chr15_5popSim_4B11_Cases_FULL.vcf.gz
```

download 2-Reference's dataset

```
mkdir anonymp/2-reference/res
cd anonymp/2-reference/res
curl -O https://lysine.univ-brest.fr/anonymp/data/chr15_5popSim_4B11_Ref.vcf.gz
```

run the automated pipeline on the 16 chunks

``` bash
cd anonymp
./pipeline.sh $(seq 16)
```
