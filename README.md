# praznik [![CRAN downloads](https://cranlogs.r-pkg.org/badges/praznik)](https://cran.rstudio.com/web/packages/praznik/index.html) 

Praznik is an R package containing a collection of information-based, feature selection filter methods.
It provides efficient implementations of several information-based filter methods, including CMIM, DISR, JMI, JMIM and mRMR, as well as corresponding feature scorers allowing for experiments.

You can install from [CRAN](https://cran.r-project.org/package=praznik), or directly from NotABug with:

```r
devtools::install_git('https://notabug.org/mbq/praznik')
```

For a bleeding edge (but working) version, install from the devel branch:

```r
devtools::install_git('https://notabug.org/mbq/praznik',branch='devel')
```

Praznik is inspired by the [FEAST library](https://github.com/Craigacp/FEAST); in fact, the earlier instance of `mbq/praznik` was an R wrapper of this software. 
This project is still available as [feast_r](https://github.com/mbq/feast_r), but is no longer maintained.

Contributions are welcome, but please make pull requests against the [devel](https://notabug.org/mbq/praznik/src/devel) branch.
