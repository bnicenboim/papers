# July 5, 2013
#    gplotlerm

library(psych)
library(knitr)
library(MASS)
library(ggplot2)
library(lme4)
#library(reshape)
library(reshape2)
library(ggplot2)
library(MASS)
#library(influence.ME)
library(scales)
library(car)
library(mixtools)
require(grid)


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

gplotlmer <- function(lmer, remeffactors,binomial=F,
                      x,group=NULL, color=NULL,linetype=NULL,shape =NULL ,facets=NULL, errorbar =F,keep=TRUE,grouping=TRUE,ran=NULL,... ) {
  
  

  vars <- c(x)
  
  #extract spaces and unwanted character from the vars
  vars <- c(vars, strsplit(paste(group,"",sep =""),split=' |:|\\+|~' )[[1]])
  vars <- c(vars, strsplit(paste(color,"",sep =""),split=' |:|\\+|~' )[[1]])
  vars <- c(vars, strsplit(paste(linetype,"",sep =""),split=' |:|\\+|~' )[[1]])
  vars <- c(vars, strsplit(paste(shape,"",sep =""),split=' |:|\\+|~' )[[1]])
  vars <- c(vars, strsplit(paste(facets,"",sep =""),split=' |:|\\+|~' )[[1]])
  vars <- unique(as.vector(vars))
  vars <- vars[!(vars %in% c("", ".")) ] #remove the empty "" and .

  
  #terms <- as.formula(paste(terms," ~." ))
  #print(vars)
  data <- model.frame(lmer)    #lmer@frame
  if(!missing( remeffactors)){
    if(!(binomial)){
    data$DV<- remef(lmer,fix=remeffactors ,keep=keep,grouping=grouping,ran=ran)
    } else {
      data$DV<- remef(lmer,fix=remeffactors ,keep=keep,grouping=grouping,ran=ran,family=binomial())
    }
      
  }else {  #not remef
    data$DV <-  data[,c(1)]
  }
  MEANS <- ddply(data,vars,summarize,M= mean(DV,na.rm=T) , SE =sd(DV,na.rm=T)/sqrt(length(DV)))
  
  
  if( !is.null(facets)   ){
    
    facet_custom <- facet_grid(as.formula(facets))
        
  } else {
    facet_custom <- NULL
  }

 
  
#  print(MEANS)
  p <- ggplot(MEANS,aes_string(x=x,y="M",group=group, color=color, linetype = linetype,shape =shape),...)+
   theme( text = (element_text(size=22)))+facet_custom 
  
  if(errorbar){
    
    p <-  p+geom_errorbar(aes(ymin = M-2*SE, ymax = M+2*SE,width=0.3))
    
  }
  return(p)
  
}

###


inv_format_show <- function() {
  function(x) format(-1000/x,digits = 2) 
}



MSE <- function(x)  c(M=mean(x), SE=sd(x)/sqrt(length(x)))

#absmax is maximum absolute value of a residual in the model, onlyone removes only the further away outlier
outliers <- function(lmer, absmax,onlyone=F ,fitted_values=F,range2keep=NULL ) {
  if(!is.null(slot(lmer,"call")$na.action) )  {
    warning("I think that na.action shouldn't be set to na.exclude, at least if the lmer is used with a subset or if there're NA")
    
  }
  #and the same data:
  data <- lmer@frame  #= slot(lmer,"frame")  
  if(is.null(range2keep)){
    range2keep <- c(-absmax,absmax)  
  }
  
  
  if(onlyone){
    if(fitted_values){
      outliersRow <- which( abs(fitted(lmer)) == max(abs( fitted(lmer))) &  abs(fitted(lmer)) > absmax)  
    }else {
    outliersRow <- which( abs(residuals(lmer)) == max(abs( residuals(lmer))) &  abs(residuals(lmer)) > absmax)
}
    
  }else {
    if(fitted_values){
      outliersRow <- which( fitted(lmer) <  range2keep[1] | fitted(lmer) >  range2keep[2] )  #residuals outside this range
      plot(fitted(lmer), resid(lmer))
      abline(v=range2keep)
    }else {
    outliersRow <- which( residuals(lmer) <  range2keep[1] | residuals(lmer) >  range2keep[2] )  #residuals outside this range
    plot(fitted(lmer), resid(lmer))
    abline(h=range2keep)
    }
   
    propOut <-  signif(length(outliersRow)/nrow(data) * 100,4)
    propOut <-  signif(length(outliersRow)/nrow(data) * 100,4)
    
    print("********************")
    print(paste( "Number of outliers:",length(outliersRow)))
    print(paste( "Proportion of outliers: ",propOut,"%"  ))
    print("*****************")
    
  }
  return(data[outliersRow,])
  
}





ci <- function(x){
  # Computes 95% confidence interval
  m <- mean(x, na.rm = TRUE)
  n <- length(x[!is.na(x)])
  
  if (length(unique(x)) > 2) {
    ## Assume t-distributed sample
    s <- sd(x, na.rm = TRUE)
    e <- qt(.975, df = n - 1) * (s / sqrt(n))
    lower <- m - e
    upper <- m + e
  } else {
    ## Proportions
    s <- sqrt((m * (1 - m)) / n)
    lower <- m - (1.96 * s)
    upper <- m + (1.96 * s)
  }
  return(data.frame(lower = lower, upper = upper))
}

changename <- function(data,oldname,newname){
   row.names(data)[which(row.names(data) == oldname)] <- newname #this is to change row name
  data
}

renaming <- function(thetable) {
  thetable <- changename(thetable,"wmc","WMC" )
  thetable <- changename(thetable,"rs","RS" )
  thetable <- changename(thetable,"cV1vsV3","VP1 vs.\\ VP3" )
  thetable <- changename(thetable,"cV2vsV1V3","VP2 vs.\\ VP1 VP3" )
  thetable <- changename(thetable,"cV2vsV1","VP2 vs.\\ VP1" )
  thetable <- changename(thetable,"cV3vsV2","VP3 vs.\\ VP2" )
  thetable <- changename(thetable,"c2V3vsV1","VP3 vs.\\ VP1" )
  thetable <- changename(thetable,"c2V2vsV3","VP2 vs.\\ VP3" )
  
  thetable <- changename(thetable,"I(wmc^2)","WMC^2" )
  thetable <- changename(thetable,"cV1vsV3:wmc","VP1 vs.\\ VP3 : WMC" )
  thetable <- changename(thetable,"cV2vsV1V3:wmc","VP2 vs.\\ VP1 VP3 : WMC" )
  thetable <- changename(thetable,"cV2vsV1:wmc","VP2 vs.\\ VP1 : WMC" )
  thetable <- changename(thetable,"cV3vsV2:wmc","VP3 vs.\\ VP2 : WMC" )
  thetable <- changename(thetable,"c2V3vsV1:wmc","VP3 vs.\\ VP1 : WMC" )
  thetable <- changename(thetable,"c2V2vsV3:wmc","VP2 vs.\\ VP3 : WMC" )
  
  thetable <- changename(thetable,"cV1vsV3:I(wmc^2)","VP1 vs.\\ VP3 : WMC^2" )
  thetable <- changename(thetable,"cV2vsV1V3:I(wmc^2)","VP2 vs.\\ VP1 VP3 : WMC^2" )
  thetable <- changename(thetable,"cV2vsV1:I(wmc^2)","VP2 vs.\\ VP1 : WMC^2" )
  thetable <- changename(thetable,"cV3vsV2:I(wmc^2)","VP3 vs.\\ VP2 : WMC^2" )
  thetable <- changename(thetable,"c2V3vsV1:I(wmc^2)","VP3 vs.\\ VP1 : WMC^2" )
  thetable <- changename(thetable,"c2V2vsV3:I(wmc^2)","VP2 vs.\\ VP3 : WMC^2" )
  
  thetable <- changename(thetable,"cV1vsV3:rs","VP1 vs.\\ VP3 : RS" )
  thetable <- changename(thetable,"cV2vsV1V3:rs","VP2 vs.\\ VP1 VP3 : RS" )
  thetable <- changename(thetable,"cV2vsV1:rs","VP2 vs.\\ VP1 : RS" )
  thetable <- changename(thetable,"cV3vsV2:rs","VP3 vs.\\ VP2 : RS" )
  thetable <- changename(thetable,"c2V3vsV1:rs","VP3 vs.\\ VP1 : RS" )
  thetable <- changename(thetable,"c2V2vsV3:rs","VP2 vs.\\ VP3 : RS" )
  
  thetable <- changename(thetable,"scale(log(FFD), scale = F)","FFD" )
  thetable <- changename(thetable,"cV1vsV3:scale(log(FFD), scale = F)","VP1 vs.\\ VP3 : FFD" )
  thetable <- changename(thetable,"cV2vsV1V3:scale(log(FFD), scale = F)","VP2 vs.\\ VP1 VP3 : FFD" )
  thetable <- changename(thetable,"cV2vsV1:scale(log(FFD), scale = F)","VP2 vs.\\ VP1 : FFD" )    
  thetable <- changename(thetable,"cV3vsV2:scale(log(FFD), scale = F)","VP3 vs.\\ VP2 : FFD" )  
  thetable <- changename(thetable,"c2V3vsV1:scale(log(FFD), scale = F)","VP3 vs.\\ VP1 : FFD" )   
  thetable <- changename(thetable,"c2V2vsV3:scale(log(FFD), scale = F)","VP2 vs.\\ VP3 : FFD" )    
  
  thetable <- changename(thetable,"cV1vsV3:wmc:scale(log(FFD), scale = F)","VP1 vs.\\ VP3 : WMC : FFD" )
  thetable <- changename(thetable,"cV2vsV1V3:wmc:scale(log(FFD), scale = F)","VP2 vs.\\ VP1 VP3 : WMC : FFD" )    
  thetable <- changename(thetable,"cV2vsV1:wmc:scale(log(FFD), scale = F)","VP2 vs.\\ VP1 : WMC : FFD" )    
  thetable <- changename(thetable,"cV3vsV2:wmc:scale(log(FFD), scale = F)","VP3 vs.\\ VP2 : WMC : FFD" )    
  thetable <- changename(thetable,"c2V3vsV1:wmc:scale(log(FFD), scale = F)","VP3 vs.\\ VP1 : WMC : FFD" )    
  thetable <- changename(thetable,"c2V2vsV3:wmc:scale(log(FFD), scale = F)","VP2 vs.\\ VP3 : WMC : FFD" )    
  
  thetable <- changename(thetable,"wmc:scale(log(FFD), scale = F)","WMC : FFD" )    
  
  thetable <- changename(thetable,"c_wordlength","word length" )
  thetable <- changename(thetable,"c_wordn","word number" )
  thetable <- changename(thetable,"wmc:c_wordlength","WMC:word length" )
  thetable <- changename(thetable,"rs:c_wordlength","RS:word length" )
}


format.glmer <- function(model,digits =1, nsmall =1){
  model <- coef(summary(model))
  
  formated <- data.frame(cbind(format.data.frame(data.frame(model)[1:3],digits =digits,scientific=F,nsmall=nsmall),format.pval(data.frame(model)[4],digits=1,nsmall=2,eps=0.001,scientific=F)))
  formated[,4] <- as.character(formated[,4])
  
  colnames(formated) <- c("Coef","\\textit{SE}","\\textit{z}","\\textit{p}")
  formated <- renaming(formated)
  formated
}

format.lmer <- function(model,digits =1, nsmall =1){
  model <- coef(summary(model))
  
  formated <- data.frame(format.data.frame(data.frame(model)[1:3],digits =digits,scientific=F,nsmall=nsmall))
  formated$signif<-ifelse(abs(as.numeric(formated$t.value))>=2," * ","   " )
  colnames(formated) <- c("Coef","\\textit{SE}","\\textit{t }","{}")
  formated <- renaming(formated)
  formated
}


tabular.table.frontiers <- function(table,preheader="",linesinrows=c() ,linesincols=c(), tightcols=c(ncol(table))){
  
  cs <-"l"
  for(c in 1:(ncol(table))){
    if(c %in%linesincols) {  #make a line after certain rows
      cs <- paste(cs,"|")
      
      
    }
    if(c %in% tightcols){ #for p-values
      cs <- paste(cs,"@{ }r@{ }")    #tight column for the p-value
      
    }else {
      cs <- paste(cs,"r")
    }
  }
  
  begin <- paste("\\begin{tabular}{",cs," }\\toprule",sep="")
  
  
  
  cat(begin)
  cat("\n")
  if(preheader!=""){
    cat(preheader)
    cat("\n")
    cat("\\midrule\\\\")
    cat("\n")
  }
  
  
  header <- paste("Predictor & ",paste(colnames(table),collapse=" & "), "\\\\","\\midrule",sep=" ")
  cat(header)
  cat("\n")
  for(i in 1:nrow(table)){
    if(i %in%linesinrows) {  #make a line after certain rows
      cat("[3pt] \\midrule \\rule{0pt}{3ex}")
      cat("\n")
    }
    
    completerow <-  paste( paste(row.names(table[i,]),  #predictor
                                 paste(table[i,],collapse=" & "),sep =" & "), #values
                           "\\\\",sep=" ")
    cat(completerow)
    cat("\n")
  }
  
  cat("\\botrule")
  cat("\n")
  
  
 # cat("\\hline")  #so that I can put two tables together
  
  cat("\\end{tabular}")
  
  cat("\n")
}



colTitlesT <- c("Coef","\\textit{SE}","\\textit{t}")
colTitlesZ  <- c("Coef","\\textit{SE}","\\textit{z}","\\textit{p}")





latexcoef <- function(model,predictor,decimals=2){
  printcoef <- ""
  values <- coef(summary(model))[ predictor,]
  if(length(values) ==4){
    coefs <- colTitlesZ
  } else {
    coefs <- colTitlesT
    
  }
  for(i in 1:length(coefs)){
    printcoef <- paste(printcoef,"\\mbox{",coefs[i]," = ", sep="")
    if(i != 4) #if it's not a p-value
    {
      printcoef <- paste(printcoef,'$',(format(round(values[i], decimals), nsmall=decimals)),"$","}",sep="")
    } else {
      printcoef <- paste(printcoef,"$",format.pval(values[i], digits=2,eps=0.0001,scientific=F ),"$","}",sep="")
      
    }
    if(i!=length(coefs)){ printcoef<- paste(printcoef,", ",sep="") }
    
  }
  
  printcoef
}


library(stringr)
ann_text <- function(ggplot,x,y,facet,text,color="black",adj=0)
{
  p<-ggplot
  fx <- as.character(p$labels$x)
  fy <- as.character(p$labels$y)
  fgroup <- as.character(p$labels$group )
  ffacet<- as.character(p$facet)[1]
  ffacet<- str_match(string=ffacet,pattern="\\((.*)\\)")[2]
  #flinetype <- as.character(p$labels$linetype )
  da<-data.frame(x, y, "", facet)
  colnames(da)<- c(fx,fy,fgroup,ffacet)
  ltype <- levels(p$data[[fgroup]])
  p + geom_text(data=da,label=text,col=color,adj=adj,aes(linetype = c("short","long") ) )  #I should fix it
}



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#
