# MTARclone
### What is MTAR:

  - an R package implementing methods for genome-wide association test of multiple traits

### How to install the original MTAR library:

  - installation within R: devtools::install_github("baolinwu/MTAR")

### Why a clone of the original package is needed:

- There is another package in CRAN named MTAR (https://cran.r-project.org/web/packages/MTAR/index.html) and when updating all R packages, the MTAR package created by Dr. Baolin Wu will be automatically expunged and changed to the package deposited on CRAN.

### Credits:

- All credits for creating and developing the MTAR package should be due to Dr. Baolin Wu. Thi is just a clone of the original repository maintained by Wanjun Gu (wagu@health.ucsd.edu).



### How to install the clone package

```R
if(!require("devtools")){
  install.packages("devtools")
  library("devtools")
}

devtools::install_github("Broccolito/MTARclone")
```



### More about the MTAR package

#### Fast and accurate genome-wide association test of multiple quantitative traits

  - Reference
    - Wu,B. and Pankow,J.S. (2018) Fast and accurate genome-wide association test of multiple continuous traits. *CMMM*, vol. 2018, Article ID 2564531, 1-9.
  - Jointly test the association of a SNP with multiple continuous traits.
    - accurate calculation of P-values.
    - very efficient and extremely scalable to genome-wide association test.
  - Sample codes

```r
  library(VGAM)
  library(MTAR)
  Z = rbinom(1000,1,0.5)
  G = rbinom(1000,2,0.25)
  ##
  X = rnorm(1000)
  e = rnorm(1000)
  Y1 = Z+X + 0.15*G + rnorm(1000)+e
  Y2 = Z+X + 0.1*G + rnorm(1000)+e
  Y = cbind(Y1,Y2)
  obj = MLM.null(Y,cbind(Z,X))
  MQTAc(obj, G)
  ##
  X1 = rnorm(1000)
  X2 = rnorm(1000)
  Y1 = Z+X1 + 0.15*G + rnorm(1000)+e
  Y2 = Z+X2 + 0.1*G + rnorm(1000)+e
  Y = cbind(Y1,Y2)
  YX = list(Y, cbind(X1,Z), cbind(X2,Z))
  objd = MDM.null(YX,pux=2)
  MQTAd(objd, G)
```



#### Multinomial ACL model for multi-trait association test

 - Reference
   - Wu,B. and Pankow,J.S. (2015) Statistical methods for association tests of multiple continuous traits in genome-wide association studies. *Annals of human genetics*, 79(4), 282-293.
 - We are testing the joint effects of multiple continuous traits with a SNP.
 - Three tests are implemented: an omnibus m-DF chi-square test; two 1-DF chi-square tests assuming common effect or common scaled effect.
 - Sample R codes

```R
library(MTAR)
n=1e3; K=4; m=3; maf=0.2
Xs = matrix(rnorm(n*K),n,K)
E = matrix(rnorm(n*m), n,m)*sqrt(0.75) + rnorm(n)*sqrt(0.25)
G = rbinom(n,2,maf)
Ys = rowMeans(Xs) + E
Y1 = Ys
Y1[,1] = Y1[,1] + 0.15*G
Y1[,2] = Y1[,2] - 0.15*G
MTA.ACL(Y1,Xs,G)
Y2 = Ys
Y2[,1] = Y2[,1] + 0.15*G
Y2[,2] = Y2[,2] + 0.15*G
MTA.ACL(Y2,Xs,G)
Y3 = Ys + 0.15*G
MTA.ACL(Y3,Xs,G)
```

