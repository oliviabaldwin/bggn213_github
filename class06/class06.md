# Class 6: Functions in R
Olivia Baldwin

First function. To create functions:
`function_name <- function(arguments){function body}`

``` r
add <- function(x,y) {
  x + y
}
```

Does the function work?

``` r
add(3,4)
```

    [1] 7

``` r
add(c(100,10,1), 1)
```

    [1] 101  11   2

Make a function “generate_dna” that makes a random nucleotide sequence
of any length.

``` r
bases <- c("A", "C", "G", "T") 
sequence <-sample(bases, size = 10, replace=TRUE)

#`replace = TRUE` allows to sample the same letter every time, i.e. it replaces the letter back into the sample pool every time
```

Above is my “snipet” that works. Now it can become a function.

``` r
generate_dna <- function(length){
  bases <- c("A", "C", "G", "T") 
  sequence <-sample(bases, size = length, replace=TRUE)
  return(sequence)
}
```

``` r
generate_dna(10)
```

     [1] "G" "A" "G" "C" "T" "C" "A" "A" "C" "G"

``` r
generate_dna(12)
```

     [1] "C" "A" "A" "G" "A" "T" "T" "A" "T" "C" "A" "T"

Lets make a protein sequence generator.

``` r
library(bio3d)
```

``` r
aa <- unique(bio3d::aa.table$aa1[1:20])
aa
```

     [1] "A" "R" "N" "D" "C" "Q" "E" "G" "H" "I" "L" "K" "M" "F" "P" "S" "T" "W" "Y"
    [20] "V"

``` r
generate_prot <- function(length){
  aa <- unique(bio3d::aa.table$aa1[1:20])
  sequence <- sample(aa, size=length, replace=TRUE)
  sequence <- paste(sequence, collapse = "")
  return(sequence)
}

#collapse puts all of the letters into one "" string
#paste will literally paste things together or paste something onto the end of each part of the vector
```

``` r
generate_prot(10)
```

    [1] "TLIRRNQCGC"

Generate random protein sequences of **length 6 to 12**. To do this use
the function `sapply()`.

``` r
prot_seqs <- sapply(6:12, generate_prot)
```

Format our sequences as fasta files.

``` r
cat(paste(">id.", 6:12, "\n", prot_seqs, sep=""), sep="\n")
```

    >id.6
    NYQIRA
    >id.7
    FFAPFRE
    >id.8
    EASPQHKD
    >id.9
    SKYCMVVWQ
    >id.10
    MAAPWCGCPS
    >id.11
    WEHELCNYQIF
    >id.12
    VKWYEHDYNHGI

``` r
# the "\n" means return to next line 
```
