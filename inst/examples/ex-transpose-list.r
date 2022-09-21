# transpose an example list (note that an index is missing)
lst <- list(1:3,3:5,c(3,7))
transpose_list(lst)
# observe that transposition is an involution
transpose_list(transpose_list(lst))
