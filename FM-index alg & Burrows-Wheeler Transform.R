#FM_index_located('GTTCCGTTCCGTTCC','GTTCC')
#Use this function to check the results.
#For other sequences and reads, please replace the sequence in the first blank and the read in the second blank.


bwt <- function (x,terminator='!') {
  #add terminator
  x <- paste0(x,terminator)
  
  #get string length
  n <- nchar(x)
  
  #define roration function
  rorate <- function (x) {
    paste0(substring(x,2), substring(x,1,1))
  }
  
  #create table to store all rotations
  tbl <- c(x,vector('character',n-1))
  
  #get all rorations
  for (i in 2:n) {
    tbl[i] <- rorate(tbl[i-1])
  }
  
  #sort by first column alphabetically
  tblorder <- order(tbl)
  tbl <- tbl[tblorder]
  sa <- c(0:(n-1))[tblorder]
  #return last column
  out <- sapply(tbl,substring,first=n,USE.NAMES=F)
  first_col <- sapply(tbl, substring, first=1, last=1, USE.NAMES=F)
  list(L=paste0(out,collapse=''),SA=sa,First=paste0(first_col, collapse=''))
 
  
}

#precaluclate skip for each letter in the alphabet set of BWT(x)
calM <- function (x) {
  x <- strsplit(x,'')[[1]]
  M <- table(x)
  cumsum(M)-M
}



Rank <- function(sequence, char, n) {
  subseq <- substring(sequence, 1, n)
  count <- sum(strsplit(subseq, NULL)[[1]] == char)
  return(count)
}



find_value <- function(char, skip_list) {
  
  return(skip_list[char])
}



FM_index_located <- function (S,R){
  
  
  
  l <- bwt(S)
  #print(l)
  n <- nchar(l)
  #print(n)
  m <- nchar(R)
  #print(m)
  first_col_value <- l$First
  
  
  #print(first_col_value)
  
  skip <- calM(first_col_value)
  #print(skip)
   
  
  Next_char <- substring(R,m,m)
  #print(paste("The next is", Next_char))
  Range_lower_bound <- find_value(Next_char,skip)
  #print(paste("The lb is", Range_lower_bound))
  Range_upper_bound <- find_value(Next_char,skip)+ Rank(first_col_value,Next_char,n)
  #print(paste("The ub is", Range_upper_bound))
  
  for(i in 2:m){
    Next_char <- substring(R,6-i,6-i)
    #print(paste("The next is", Next_char))
    rank_a <- Rank(l,Next_char,Range_lower_bound)
    
    Range_lower_bound=rank_a+ find_value(Next_char,skip)
    #print(paste("The lb is", Range_lower_bound))
    
    rank_b <- Rank(l,Next_char,Range_upper_bound)
    Range_upper_bound=rank_b+ find_value(Next_char,skip)
    #print(paste("The ub is", Range_upper_bound))
  }
  #print(paste("The lb is", Range_lower_bound))
  #print(paste("The ub is", Range_upper_bound))
 
    result <- seq(Range_lower_bound,Range_upper_bound-1)
    #print(result)
    #print(paste("The result is", result))
    result_number <- length(result)0
    print(paste("The number of matches is", result_number))
    for(t in 1:result_number){
      
      position <- l$SA[result[t]+1]
      print(paste("The match is at offset", position))
      
      
    }
    
  
  
}