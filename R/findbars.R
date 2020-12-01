expandDoubleVerts <- function (term) 
{
  expandDoubleVert <- function(term) {
    frml <- formula(substitute(~x, list(x = term[[2]])))
    newtrms <- paste0("0+", attr(terms(frml), "term.labels"))
    if (attr(terms(frml), "intercept") != 0) 
      newtrms <- c("1", newtrms)
    as.formula(paste("~(", paste(vapply(newtrms, function(trm) paste0(trm, 
                                                                      "|", deparse(term[[3]])), ""), collapse = ")+("), 
                     ")"))[[2]]
  }
  if (!is.name(term) && is.language(term)) {
    if (term[[1]] == as.name("(")) {
      term[[2]] <- expandDoubleVerts(term[[2]])
    }
    stopifnot(is.call(term))
    if (term[[1]] == as.name("||")) 
      return(expandDoubleVert(term))
    term[[2]] <- expandDoubleVerts(term[[2]])
    if (length(term) != 2) {
      if (length(term) == 3) 
        term[[3]] <- expandDoubleVerts(term[[3]])
    }
  }
  term
}



findbars <- function (term) 
{
  fb <- function(term) {
    if (is.name(term) || !is.language(term)) 
      return(NULL)
    if (term[[1]] == as.name("(")) 
      return(fb(term[[2]]))
    stopifnot(is.call(term))
    if (term[[1]] == as.name("|")) 
      return(term)
    if (length(term) == 2) 
      return(fb(term[[2]]))
    c(fb(term[[2]]), fb(term[[3]]))
  }
  expandSlash <- function(bb) {
    makeInteraction <- function(x) {
      if (length(x) < 2) 
        return(x)
      trm1 <- makeInteraction(x[[1]])
      trm11 <- if (is.list(trm1)) 
        trm1[[1]]
      else trm1
      list(substitute(foo:bar, list(foo = x[[2]], bar = trm11)), 
           trm1)
    }
    slashTerms <- function(x) {
      if (!("/" %in% all.names(x))) 
        return(x)
      if (x[[1]] != as.name("/")) 
        stop("unparseable formula for grouping factor", 
             call. = FALSE)
      list(slashTerms(x[[2]]), slashTerms(x[[3]]))
    }
    if (!is.list(bb)) 
      expandSlash(list(bb))
    else unlist(lapply(bb, function(x) {
      if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]]))) 
        lapply(unlist(makeInteraction(trms)), function(trm) substitute(foo | 
                                                                         bar, list(foo = x[[2]], bar = trm)))
      else x
    }))
  }
  modterm <- expandDoubleVerts(if (is(term, "formula")) 
    term[[length(term)]]
    else term)
  expandSlash(fb(modterm))
}
