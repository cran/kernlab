
setGeneric("as.kernelMatrix",function(x, ...) standardGeneric("as.kernelMatrix"))
setMethod("as.kernelMatrix", signature(x = "matrix"),
function(x, center = FALSE...)
{
  return(new("kernelMatrix",.Data = x))
})
