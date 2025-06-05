# class06
Hailey Heirigs (PID: A16962278)

## Quarto

Quarto enables you to weave together content and executable code into a
finished document. To learn more about Quarto see <https://quarto.org>.

## Running Code

When you click the **Render** button a document will be generated that
includes both content and the output of embedded code. You can embed
code like this:

``` r
1 + 1
```

    [1] 2

``` r
add <- function(x, y){
  x+y
}
```

I can just use this function like any other function as long as R knows
about it (i.e. run the code chunk)

``` r
add(1, 100)
```

    [1] 101

## 1. R functions

Every R function has 3 things. - name (we get to pick this) - input
arguments - the body

``` r
add <- function(x, y=10, z=0){
  x+y+z
}
```

``` r
add(1)
```

    [1] 11

Functions can have “required” input arguments and “optional” input
arguments. The optional arguments are defined with an equals default
value (`y=10`) in the function definition.

``` r
add(x=1, y=100, z=10)
```

    [1] 111

> Q. Write a function to return a DNA sequence of a user specified
> length? Call it `generate_dna()`

The `sample()` function can help here

``` r
#generate_dna <- function(size=5) { }

students <- c("jeff", "jeremy", "peter")

sample(students, size = 5, replace=TRUE)
```

    [1] "jeremy" "jeremy" "peter"  "peter"  "peter" 

## 2. Generate DNA sequences

Now work with `bases` rather than `students`

``` r
bases <- c("A", "C", "G", "T")

sample(bases, size = 10, replace = TRUE)
```

     [1] "T" "T" "C" "C" "G" "A" "C" "T" "C" "C"

Now I have a working snippet of code & the body of my first function is
here.

``` r
generate_dna <- function(size=5) {
  bases <- c("A", "C", "G", "T")
  sample(bases, size=size, replace=TRUE)
}
```

``` r
generate_dna(100)
```

      [1] "C" "T" "G" "C" "C" "T" "T" "A" "T" "G" "C" "G" "T" "T" "G" "T" "A" "C"
     [19] "G" "G" "T" "T" "C" "A" "A" "G" "T" "G" "G" "G" "A" "A" "C" "A" "T" "G"
     [37] "C" "T" "A" "C" "A" "T" "C" "T" "T" "G" "A" "T" "G" "G" "A" "A" "C" "G"
     [55] "G" "A" "T" "C" "G" "T" "A" "T" "A" "C" "A" "G" "G" "G" "A" "C" "G" "G"
     [73] "G" "G" "T" "G" "G" "G" "G" "A" "T" "C" "C" "A" "C" "C" "T" "T" "C" "G"
     [91] "G" "A" "G" "A" "T" "T" "C" "G" "C" "G"

``` r
generate_dna()
```

    [1] "C" "A" "G" "C" "C"

I want the ability to return a sequence like “AGTACCTG” i.e. a one
element vector where the bases are all together.

``` r
generate_dna <- function(size=5, together=TRUE) {
  bases <- c("A", "C", "G", "T")
  sequence <- sample(bases, size=size, replace=TRUE)
  if(together) {
    sequence <- paste(sequence, collapse = "")
  }
  return(sequence)
}
```

``` r
generate_dna()
```

    [1] "TATGC"

``` r
generate_dna(together=FALSE)
```

    [1] "T" "T" "C" "C" "T"

## 3. Generate Protein function

> Q. Write a protein sequence generating function that will return
> sequences of a user specified length?

We can get the set of 20 natural amino-acids from the **bio3d** package.

``` r
aa <- bio3d::aa.table$aa1[1:20]
```

``` r
aa
```

     [1] "A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K" "M" "F" "P" "S" "T" "W" "Y"
    [20] "V"

and use this in our function

``` r
generate_protein <- function(size=6, together=TRUE) {
  ## Get the 20 amino-acids as a vector
  aa <- bio3d::aa.table$aa1[1:20]
  sequence <- sample(aa, size, replace=TRUE)
  
  ## Optionally return a single element string 
  if(together){
    sequence <- paste(sequence, collapse = "")
  }
  return(sequence)
}
```

> Q. Generate random protein sequences of length 6 to 12 amino acids.

``` r
generate_protein(7)
```

    [1] "RRQPCLM"

``` r
generate_protein(8)
```

    [1] "NQFSHIGS"

``` r
generate_protein(9)
```

    [1] "RGFHWPVPC"

We can fix this inability to generate multiple sequences by either
editing and adding to the function body code (e.g. a for loop) or by
using the R **apply** family of utility functions.

``` r
sapply(6:12, generate_protein)
```

    [1] "THCPCN"       "NYLCIWY"      "SDHTWHVP"     "EWSDELVWL"    "VREEQYVFQY"  
    [6] "SNPGQPASKNQ"  "SDGNVDRKMPDS"

> Q. Determine if these sequences can be found in nature or are they
> unique? Why or why not?

It would be cool and useful if I could get FASTA format output

``` r
ans <- sapply(6:12, generate_protein)
ans
```

    [1] "PCPCCE"       "ETKYCRT"      "TKNQLRQR"     "WVNAKLICC"    "IRQQKTIYML"  
    [6] "FNPGWSIHLSV"  "LYCCVTGRIPTF"

``` r
cat(ans, sep="\n")
```

    PCPCCE
    ETKYCRT
    TKNQLRQR
    WVNAKLICC
    IRQQKTIYML
    FNPGWSIHLSV
    LYCCVTGRIPTF

I want this to look like

    >ID.6
    QCRAWKLHID
    >ID.7
    TYNVGCMDFNT
    >ID.8
    EVLPAFWHLWIV

The functions `paste()` and cat()\` can help us here…

``` r
ans2 <-paste(">ID.", 7:12, "\n",ans, sep="")
cat(ans2, sep = "\n")
```

    >ID.7
    PCPCCE
    >ID.8
    ETKYCRT
    >ID.9
    TKNQLRQR
    >ID.10
    WVNAKLICC
    >ID.11
    IRQQKTIYML
    >ID.12
    FNPGWSIHLSV
    >ID.7
    LYCCVTGRIPTF

> Q. Determine if these sequences can be found in nature or are they
> unique? Why or why not?

I BLASTp searched my FASTA format sequences against NR and found that
lengths 6, 7, 8, are not unique and can be found in the databases with
100% coverage and 100% identity.

Random sequences of length 9 and above are unique and can’t be found in
the databases.

You can add options to executable code like this

    [1] 4

The `echo: false` option disables the printing of code (only output is
displayed).

``` r
1+1
```

    [1] 2
