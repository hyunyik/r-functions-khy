cap_first_letter <- function(strings, exception.list = NULL) {
  if(!require(stringr)) install.packages("stringr")
  
  first_word <- strsplit(strings, " ")[[1]][1]
  #print(first_word)
  if(first_word %in% exception.list) {
    return(strings)
  } else {
    return(str_to_sentence(strings))
  }
}

insert_new_line <- function(strings, max.num = 30) {
  return(sapply(strings, function(x) {
    paste(strwrap(x, max.num), collapse = "\n")
  }))
}

process_go <- function(strings, max.num = 30, cap.exp = NULL, new_line = T) {
#  r <- sapply(strings, function(x) {
#    ifelse(new_line, 
#           paste(strwrap(cap_first_letter(x, cap.exp), max.num), collapse = "\n"),
#           strtrim(cap_first_letter(x, cap.exp), max.num))
#  })

  r <- NA
  
  if(new_line) {
    r <- sapply(strings, function(x) { paste(strwrap(cap_first_letter(x, cap.exp), max.num), collapse = "\n") })
  } else {
    r <- sapply(strings, function(x) {
      print(paste0(x, " ",nchar(x)))
      ifelse(nchar(x)>max.num, paste0(strtrim(cap_first_letter(x, cap.exp), max.num), " ..."), cap_first_letter(x, cap.exp))
      })
  }
  
  return(r)
}
