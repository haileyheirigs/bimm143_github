# Class 7: Machine Learning 1
Hailey Heirigs (PID: A16962278)

- [Clustering](#clustering)
  - [K-means](#k-means)
  - [Heirarchical Clustering](#heirarchical-clustering)
- [Principal Component Analysis
  (PCA)](#principal-component-analysis-pca)
  - [Data Import](#data-import)
  - [PCA to the rescue](#pca-to-the-rescue)

Today we will explore unsupervised machine learning methods starting
with clustering and dimensionality reduction.

## Clustering

To start let’s make up some data to cluster where we know what the
answer should be. The `rnorm()` function will help us here.

``` r
hist( rnorm(10000, mean=3) )
```

![](class07_files/figure-commonmark/unnamed-chunk-1-1.png)

Return 30 numbers cenetered on -3 and +3

``` r
tmp <- c( rnorm(30, mean=-3),
        rnorm(30, mean=+3) )

x <- cbind(x=tmp, y=rev(tmp))

x
```

                   x          y
     [1,] -2.5494679  2.1415777
     [2,] -4.3639633  4.4896195
     [3,] -2.1829470  2.5401347
     [4,] -2.6038048  2.7126408
     [5,] -3.3818276  3.0461225
     [6,] -1.5809164  1.6558613
     [7,] -2.2520406  1.2001408
     [8,] -2.8634059  2.2120638
     [9,] -2.1306104  3.6015440
    [10,] -2.9444684  1.8654866
    [11,] -3.4391723  2.5846914
    [12,] -2.2429070  2.3445938
    [13,] -3.1811296  2.7563866
    [14,] -3.2202585  3.2276449
    [15,] -1.0978797  4.1807797
    [16,] -2.2130921  3.1815989
    [17,] -3.8698691  3.2223958
    [18,] -0.4760691  2.2033150
    [19,] -1.2317802  2.8917213
    [20,] -2.7048102  3.8934054
    [21,] -1.4907738  3.0675433
    [22,] -2.7515811  2.2287164
    [23,] -2.5588170  1.5640124
    [24,] -6.1679689  3.0254873
    [25,] -2.0700350  4.2675982
    [26,] -2.6601255  3.4946760
    [27,] -2.9307086  3.7646258
    [28,] -3.7748418  4.2650705
    [29,] -2.5833935  3.0980992
    [30,] -2.4149062  2.6033542
    [31,]  2.6033542 -2.4149062
    [32,]  3.0980992 -2.5833935
    [33,]  4.2650705 -3.7748418
    [34,]  3.7646258 -2.9307086
    [35,]  3.4946760 -2.6601255
    [36,]  4.2675982 -2.0700350
    [37,]  3.0254873 -6.1679689
    [38,]  1.5640124 -2.5588170
    [39,]  2.2287164 -2.7515811
    [40,]  3.0675433 -1.4907738
    [41,]  3.8934054 -2.7048102
    [42,]  2.8917213 -1.2317802
    [43,]  2.2033150 -0.4760691
    [44,]  3.2223958 -3.8698691
    [45,]  3.1815989 -2.2130921
    [46,]  4.1807797 -1.0978797
    [47,]  3.2276449 -3.2202585
    [48,]  2.7563866 -3.1811296
    [49,]  2.3445938 -2.2429070
    [50,]  2.5846914 -3.4391723
    [51,]  1.8654866 -2.9444684
    [52,]  3.6015440 -2.1306104
    [53,]  2.2120638 -2.8634059
    [54,]  1.2001408 -2.2520406
    [55,]  1.6558613 -1.5809164
    [56,]  3.0461225 -3.3818276
    [57,]  2.7126408 -2.6038048
    [58,]  2.5401347 -2.1829470
    [59,]  4.4896195 -4.3639633
    [60,]  2.1415777 -2.5494679

Make a plot of `x`

``` r
plot(x)
```

![](class07_files/figure-commonmark/unnamed-chunk-3-1.png)

### K-means

The main function in “base” R for K-means clustering is called
`kmeans()`:

``` r
km <- kmeans(x, centers = 2)
km 
```

    K-means clustering with 2 clusters of sizes 30, 30

    Cluster means:
              x         y
    1  2.911030 -2.664452
    2 -2.664452  2.911030

    Clustering vector:
     [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1
    [39] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

    Within cluster sum of squares by cluster:
    [1] 53.93439 53.93439
     (between_SS / total_SS =  89.6 %)

    Available components:

    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

The `kmeans()` function return a “list” with 9 components. You can see
the named components of any list with the `attributes` function.

``` r
attributes(km)
```

    $names
    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

    $class
    [1] "kmeans"

> Q. How many points are in each cluster?

``` r
km$size
```

    [1] 30 30

> Q. Cluster assignment/membership vector?

``` r
km$cluster
```

     [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1
    [39] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

> Q. Cluster centers?

``` r
km$centers
```

              x         y
    1  2.911030 -2.664452
    2 -2.664452  2.911030

> Q. Make a plot of our `kmeans()` results showing cluster assignment
> using different colors for each cluster/group of points and cluster
> centers in blue.

``` r
plot(x, col=km$cluster )
points(km$centers, col="blue", pch=15, cex=2)
```

![](class07_files/figure-commonmark/unnamed-chunk-9-1.png)

> Q. Run `kmeans()` again on `x` and this time cluster into 4
> groups/clusters and plot the same result figure as above.

``` r
km4 <- kmeans(x, centers = 4)
km4 
```

    K-means clustering with 4 clusters of sizes 8, 10, 30, 12

    Cluster means:
              x         y
    1  3.327177 -3.924879
    2  3.544159 -2.111321
    3 -2.664452  2.911030
    4  2.105991 -2.285111

    Clustering vector:
     [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 2 1 2 2 2 1 4
    [39] 4 2 2 2 4 1 2 2 1 1 4 1 4 2 4 4 4 1 4 4 1 4

    Within cluster sum of squares by cluster:
    [1] 10.12915  5.83406 53.93439  7.30444
     (between_SS / total_SS =  92.6 %)

    Available components:

    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

``` r
attributes(km4)
```

    $names
    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

    $class
    [1] "kmeans"

``` r
km4$size 
```

    [1]  8 10 30 12

``` r
km4$cluster
```

     [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 2 1 2 2 2 1 4
    [39] 4 2 2 2 4 1 2 2 1 1 4 1 4 2 4 4 4 1 4 4 1 4

``` r
km4$centers 
```

              x         y
    1  3.327177 -3.924879
    2  3.544159 -2.111321
    3 -2.664452  2.911030
    4  2.105991 -2.285111

``` r
plot(x, col=km4$cluster )
points(km4$centers, col="blue", pch=15, cex=2)
```

![](class07_files/figure-commonmark/unnamed-chunk-10-1.png)

> **key-point**: K-means clustering is super popular but can be
> mis-used. One big limitation is that it can impose a clustering
> pattern on your data even if clear natural grouping doesn’t exist -
> i.e. it does what you tell it to do in terms of `centers`.

### Heirarchical Clustering

The main function in “base” R for Heirarchical clustering is called
`hclust()`. You must give “distance matrix” as input. We can get this
from the `dist()` function in R.

You can’t just pass our dataset as is into `hclust()`

``` r
d <- dist(x)
hc <- hclust(d)
hc
```


    Call:
    hclust(d = d)

    Cluster method   : complete 
    Distance         : euclidean 
    Number of objects: 60 

The results of `hclust()` don’t have a useful `print()` method but do
have a special `plot()` method.

``` r
plot(hc)
```

![](class07_files/figure-commonmark/unnamed-chunk-12-1.png)

To get our main cluster assignment (membership vector) we need to “cut”
the tree at the big goal posts…

``` r
cutree(hc, h=8)
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
grps <- cutree(hc, h=8)
grps
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
table(grps)
```

    grps
     1  2 
    30 30 

``` r
plot(x, col=grps)
```

![](class07_files/figure-commonmark/unnamed-chunk-16-1.png)

Hierarchical Clustering is distinct in that the dendrogram (tree figure)
can reveal the potential grouping in your data (unlike K-means).

## Principal Component Analysis (PCA)

PCA is a common and highly useful dimensionality reduction technique
used in many fields - particularly bioinformatics.

Here we will analyze some data from the UK on food consumption.

### Data Import

``` r
url <- "https://tinyurl.com/UK-foods"
x <- read.csv(url)

head(x)
```

                   X England Wales Scotland N.Ireland
    1         Cheese     105   103      103        66
    2  Carcass_meat      245   227      242       267
    3    Other_meat      685   803      750       586
    4           Fish     147   160      122        93
    5 Fats_and_oils      193   235      184       209
    6         Sugars     156   175      147       139

``` r
rownames(x) <- x[,1]
x <- x[,-1]
head(x)
```

                   England Wales Scotland N.Ireland
    Cheese             105   103      103        66
    Carcass_meat       245   227      242       267
    Other_meat         685   803      750       586
    Fish               147   160      122        93
    Fats_and_oils      193   235      184       209
    Sugars             156   175      147       139

``` r
x <- read.csv(url, row.names = 1)
head(x)
```

                   England Wales Scotland N.Ireland
    Cheese             105   103      103        66
    Carcass_meat       245   227      242       267
    Other_meat         685   803      750       586
    Fish               147   160      122        93
    Fats_and_oils      193   235      184       209
    Sugars             156   175      147       139

``` r
barplot(as.matrix(x), beside=F, col=rainbow(nrow(x)))
```

![](class07_files/figure-commonmark/unnamed-chunk-20-1.png)

One conventional plot that can be useful is called a “pairs” plot.

``` r
pairs(x, col=rainbow(nrow(x)), pch=16)
```

![](class07_files/figure-commonmark/unnamed-chunk-21-1.png)

### PCA to the rescue

The main function in base R for PCA is called `prcomp()`.

``` r
pca <- prcomp( t(x) )
summary(pca)
```

    Importance of components:
                                PC1      PC2      PC3       PC4
    Standard deviation     324.1502 212.7478 73.87622 2.921e-14
    Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

The `prcomp()` function returns a list object of our results with five
attributes/components

``` r
attributes(pca)
```

    $names
    [1] "sdev"     "rotation" "center"   "scale"    "x"       

    $class
    [1] "prcomp"

The two main “results” in here are `pca$x` and `pca$rotation`. The first
of these (`pca$x`) contains the scores of the data on the new PC axis -
we use these to make our “PCA plot”.

``` r
pca$x
```

                     PC1         PC2        PC3           PC4
    England   -144.99315   -2.532999 105.768945 -9.152022e-15
    Wales     -240.52915 -224.646925 -56.475555  5.560040e-13
    Scotland   -91.86934  286.081786 -44.415495 -6.638419e-13
    N.Ireland  477.39164  -58.901862  -4.877895  1.329771e-13

``` r
library(ggplot2)
library(ggrepel)

# Make a plot of pca$x with PC1 vs PC2 
ggplot(pca$x) + 
  aes(PC1, PC2, label=rownames(pca$x)) +
  geom_point() + 
  geom_text_repel() 
```

![](class07_files/figure-commonmark/unnamed-chunk-25-1.png)

The second major result is contained in the `pca$rotation` object or
component. Let’s plot this to see what PCA is picking up…

``` r
ggplot(pca$rotation) + 
  aes(PC1, rownames(pca$rotation)) +
  geom_col()
```

![](class07_files/figure-commonmark/unnamed-chunk-26-1.png)

``` r
ggplot(pca$rotation) + 
  aes(PC2, rownames(pca$rotation)) +
  geom_col()
```

![](class07_files/figure-commonmark/unnamed-chunk-27-1.png)

These plots show comparisons between the different countries and their
food eating habits. It compares the different kinds of foods and helps
us to understand their individual health habits as well. We saw that
Ireland was an outlier when making a scatter plot from ggplot2, as well
as plotted and compared the values from PC1 and PC2 as dimensionality
reduction techniques to see where the differences deviate and are more
simlar. The generated bar plots from using geom_col helped to illustrate
what stuck out strongly for comparison as well.
