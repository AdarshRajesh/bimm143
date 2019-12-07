Packages
================
Adarsh
10/22/2019

## R Functions Revisited

Source my functions from last day

``` r
source("http://tinyurl.com/rescale-R")
```

``` r
x <- c(1,2,NA,3,NA)
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
both_na <- function(x,y) {
  sum(is.na(x) & is.na(y))
}
```

``` r
x <- c(1,NA,NA)
y1 <- c(1,NA,NA)
y2 <-c(1,NA,NA,NA,NA,NA,NA)
which(is.na(y2))
```

    ## [1] 2 3 4 5 6 7

``` r
both_na(x,y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 4

``` r
both_na <- function(x,y) {
  if(length(x)!=length(y)){
    stop("Inputs x and y should be same length")
  }
  sum(is.na(x) & is.na(y))
}
```

``` r
x <- c(1,NA,NA)
y1 <- c(1,NA,NA)
y2 <-c(1,NA,NA,NA,NA,NA,NA)
#both_na(x,y2)
```

``` r
grades <- function(x){
  if(sum(is.na(x))==0){
    min  = which.min(x)
    total = sum(x[-min])
  }
  else{
    total = sum(x,na.rm = TRUE)
  }
  grade = total/(length(x)-1)
  return(grade)
}
```

``` r
x<- c(100,100,100,NA,NA)
grades(x)
```

    ## [1] 75

``` r
url<- "https://tinyurl.com/gradeinput"
hw <- read.csv(url,row.names = 1)
```

``` r
apply(hw,1,grades)
```

    ##  student-1  student-2  student-3  student-4  student-5  student-6 
    ##      91.75      82.50      84.25      84.25      88.25      89.00 
    ##  student-7  student-8  student-9 student-10 student-11 student-12 
    ##      94.00      93.75      87.75      79.00      86.00      91.75 
    ## student-13 student-14 student-15 student-16 student-17 student-18 
    ##      92.25      87.75      78.75      89.50      88.00      94.50 
    ## student-19 student-20 
    ##      82.75      82.75
