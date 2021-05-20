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

process_go2 <- function(strings, new_line.max = 30, cut.max = 30, cap.exp = NULL, type = "new_line", verbose = F) {
  
  r <- NA
  
  if(type == "new_line") {
    if(verbose) print(paste0("type = ", type))
    r <- sapply(strings, function(x) { paste(strwrap(cap_first_letter(x, cap.exp), new_line.max), collapse = "\n") })
  } else if(type =="cut") {
    if(verbose) print(paste0("type = ", type))
    r <- sapply(strings, function(x) {
      #print(paste0(x, " ",nchar(x)))
      ifelse(nchar(x)>cut.max, paste0(strwrap(cap_first_letter(x, cap.exp), cut.max)[1], " ..."), cap_first_letter(x, cap.exp))
    })
  } else if(type =="both") {
    if(verbose) print(paste0("type = ", type))
    if(new_line.max>cut.max) {
      print("process_go2: new_line.max can't be bigger than cut.max.")
      return(FALSE)
    }
    r <- sapply(strings, function(x) {
      if(verbose) print(paste0("cut.max = ", cut.max, "| length = ", nchar(x), "| string = ", x))
      if(nchar(x)>cut.max) {
        xx <- paste0(strwrap(cap_first_letter(x, cap.exp), cut.max)[1], " ...")
        paste(strwrap(xx, new_line.max), collapse = "\n")
      } else {
        paste(strwrap(cap_first_letter(x, cap.exp), new_line.max), collapse = "\n")
      }
    })
    
  }
  
  for(u in unique(r[duplicated(r)])) {
    i <- 1
    for(idx in which(r == u)) {
      if(verbose) print(paste0("Process for duplicates: ", u, " (", i, ")"))
      r[idx] <- paste0(u, " (", i, ")")
      i <- i + 1
    }
  }
  
  return(r)
}