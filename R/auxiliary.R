#Auxiliary functions

intermaker <- function(f, blank) {
  # f is the creator function like 'Fiksel'
  class(f) <- c("intermaker", class(f))
  # blank is the prototype interaction object: extract some fields
  desired <- c("creator", "name", "par", "parnames", "pardesc")
  avail <- desired[desired %in% names(blank)]
  attr(f, "b") <- blank[avail]
  return(f)
}