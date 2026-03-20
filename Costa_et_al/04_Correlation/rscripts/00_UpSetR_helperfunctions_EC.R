#Source: https://github.com/hms-dbmi/UpSetR/issues/85 See docmanny's contribution

#Helper Functions for UpSetR
fromList2 <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

get_intersect_members <- function (x, ...){
  require(dplyr)
  require(tibble)
  x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
  n <- names(x)
  x %>% rownames_to_column() -> x
  l <- c(...)
  a <- intersect(names(x), l)
  ar <- vector('list',length(n)+1)
  ar[[1]] <- x
  i=2
  for (item in n) {
    if (item %in% a){
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  do.call(filter_, ar) %>% column_to_rownames() -> x
  return(x)
}




#####Testing funcitons work
#Example list:
lst <- list(a=c("CARD11_0","EZH2_0","HOXD11_0","FGFR1_0","FGFR1_1"), b=c("EZH2_0","EZH2_0","HOXD11_0","FGFR1_0","FGFR1_0"))

# Binary table with colnames:
fromList2(lst)


# get_intersect_members() takes as arguments a dataframe that has been formatted as a binary 
#     table, such as movies from the UpSetR vignette; as well as a series of strings with the names of 
#     columns you wish to test membership for.

movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")


# Find all movies that are exclusively dramas:
get_intersect_members(movies, 'Drama')

get_intersect_members(fromList(lst), 'a')  
