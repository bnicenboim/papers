unscale <- function(scaledvalue,mean,sd){(scaledvalue* sd) + mean}

nicedecimals <- function(number, decimals){  #to have nice ending decimals for graphs
  round(number*10^decimals/5,0)*5/10^decimals
}


#ggplot either the mean of the DV as specified in the lmer or the mean of the remef function (keep=TRUE,grouping=TRUE by default)

numbers2words <- function(x){
  
  helper <- function(x){
    
    digits <- rev(strsplit(as.character(x), "")[[1]])
    nDigits <- length(digits)
    if (nDigits == 1) as.vector(ones[digits])
    else if (nDigits == 2)
      if (x <= 19) as.vector(teens[digits[1]])
    else trim(paste(tens[digits[2]],
                    Recall(as.numeric(digits[1])),sep="-"))
    else if (nDigits == 3) trim(paste(ones[digits[3]], "hundred", 
                                      Recall(makeNumber(digits[2:1]))))
    else {
      nSuffix <- ((nDigits + 2) %/% 3) - 1
      if (nSuffix > length(suffixes)) stop(paste(x, "is too large!"))
      trim(paste(Recall(makeNumber(digits[
        nDigits:(3*nSuffix + 1)])),
                 suffixes[nSuffix],  
                 Recall(makeNumber(digits[(3*nSuffix):1])) ))
    }
  }
  trim <- function(text){
    gsub("^\ ", "", gsub("\ *$", "", text))
  }      
  
  
  makeNumber <- function(...) as.numeric(paste(..., collapse=""))     
  opts <- options(scipen=100) 
  on.exit(options(opts)) 
  ones <- c("", "one", "two", "three", "four", "five", "six", "seven",
            "eight", "nine") 
  names(ones) <- 0:9 
  teens <- c("ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
             
             "sixteen", " seventeen", "eighteen", "nineteen")     
  names(teens) <- 0:9 
  tens <- c("twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty",
            
            "ninety") 
  names(tens) <- 2:9 
  x <- round(x) 
  suffixes <- c("thousand", "million", "billion", "trillion")    
  
  if (length(x) > 1) return(sapply(x, helper))    
  helper(x) 
  
    
}


proper <- function(text){
  initial <- toupper(substr(text,1,1)  )
  rest  <- substr(text,2,nchar(text))
  return(paste(initial, rest,sep=""))
}
