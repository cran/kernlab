## S4 object definitions and assigment/accessor functions for the slots.
## should be done using virtual classes but current S4/Namespace conflict prevents it
##
## created  10.09.03 alexandros
## updated  15.09.04 

##Support Vector Machine object 
setClass("ksvm", representation(type = "character",
                                param = "list",
                                kernelf = "function",
                                kpar = "list",
                                kcall = "ANY",
                                scaling = "ANY",
                                kterms = "ANY",
##                                margin = "vector",
                                xmatrix = "matrix",
                                ymatrix = "ANY",
                                fit = "ANY",
                                lev = "vector",
                                nclass = "numeric",
                                prior = "list",
                                alpha = "ANY",
                                coeff = "ANY",
                                alphaindex = "ANY",
                                b = "numeric",
                                prob.model = "list",
                                SVindex = "vector",
                                nSV = "numeric",
                                error = "numeric",
                                cross = "numeric",
                                n.action= "ANY"))

## Create accessors for all class slots
## can this be done in a neater way ?

if(!isGeneric("type")){
  if (is.function("type"))
    fun <- type
  else fun <- function(object) standardGeneric("type")
  setGeneric("type", fun)
}
setMethod("type", "ksvm", function(object) object@type)
setGeneric("type<-", function(x, value) standardGeneric("type<-"))
setReplaceMethod("type", "ksvm", function(x, value) {
  x@type <- value
  x
})

if(!isGeneric("margin")){
  if (is.function("margin"))
    fun <- margin
  else fun <- function(object) standardGeneric("margin")
  setGeneric("margin", fun)
}
setMethod("margin", "ksvm", function(object) object@margin)
setGeneric("margin<-", function(x, value) standardGeneric("margin<-"))
setReplaceMethod("margin", "ksvm", function(x, value) {
  x@margin <- value
  x
})

if(!isGeneric("param")){
  if (is.function("param"))
    fun <- param
  else fun <- function(object) standardGeneric("param")
  setGeneric("param", fun)
}
setMethod("param", "ksvm", function(object) object@param)
setGeneric("param<-", function(x, value) standardGeneric("param<-"))
setReplaceMethod("param", "ksvm", function(x, value) {
  x@param <- value
  x
})

if(!isGeneric("kernelf")){
  if (is.function("kernelf"))
    fun <- kernelf
  else fun <- function(object) standardGeneric("kernelf")
  setGeneric("kernelf", fun)
}
setMethod("kernelf", "ksvm", function(object) object@kernelf)
setGeneric("kernelf<-", function(x, value) standardGeneric("kernelf<-"))
setReplaceMethod("kernelf", "ksvm", function(x, value) {
  x@kernelf <- value
  x
})

if(!isGeneric("kpar")){
  if (is.function("kpar"))
    fun <- kpar
  else fun <- function(object) standardGeneric("kpar")
  setGeneric("kpar", fun)
}
setMethod("kpar", "ksvm", function(object) object@kpar)
setGeneric("kpar<-", function(x, value) standardGeneric("kpar<-"))
setReplaceMethod("kpar", "ksvm", function(x, value) {
  x@kpar <- value
  x
})


if(!isGeneric("kcall")){
  if (is.function("kcall"))
    fun <- kcall
  else fun <- function(object) standardGeneric("kcall")
  setGeneric("kcall", fun)
}
setMethod("kcall", "ksvm", function(object) object@kcall)
setGeneric("kcall<-", function(x, value) standardGeneric("kcall<-"))
setReplaceMethod("kcall", "ksvm", function(x, value) {
  x@kcall <- value
  x
})

if(!isGeneric("scaling")){
  if (is.function("scaling"))
    fun <- scaling
  else fun <- function(object) standardGeneric("scaling")
  setGeneric("scaling", fun)
}
setMethod("scaling", "ksvm", function(object) object@scaling)
setGeneric("scaling<-", function(x, value) standardGeneric("scaling<-"))
setReplaceMethod("scaling", "ksvm", function(x, value) {
  x@scaling<- value
  x
})



if(!isGeneric("kterms")){
  if (is.function("kterms"))
    fun <- kterms
  else fun <- function(object) standardGeneric("kterms")
  setGeneric("kterms", fun)
}
setMethod("kterms", "ksvm", function(object) object@kterms)
setGeneric("kterms<-", function(x, value) standardGeneric("kterms<-"))
setReplaceMethod("kterms", "ksvm", function(x, value) {
  x@kterms <- value
  x
})

if(!isGeneric("xmatrix")){
  if (is.function("xmatrix"))
    fun <- xmatrix
  else fun <- function(object) standardGeneric("xmatrix")
  setGeneric("xmatrix", fun)
}
setMethod("xmatrix", "ksvm", function(object) object@xmatrix)
setGeneric("xmatrix<-", function(x, value) standardGeneric("xmatrix<-"))
setReplaceMethod("xmatrix", "ksvm", function(x, value) {
  x@xmatrix <- value
  x
})

if(!isGeneric("ymatrix")){
  if (is.function("ymatrix"))
    fun <- ymatrix
  else fun <- function(object) standardGeneric("ymatrix")
  setGeneric("ymatrix", fun)
}
setMethod("ymatrix", "ksvm", function(object) object@ymatrix)
setGeneric("ymatrix<-", function(x, value) standardGeneric("ymatrix<-"))
setReplaceMethod("ymatrix", "ksvm", function(x, value) {
  x@ymatrix <- value
  x
})

if(!isGeneric("fit")){
  if (is.function("fit"))
    fun <- fit
  else fun <- function(object) standardGeneric("fit")
  setGeneric("fit", fun)
}
setMethod("fit", "ksvm", function(object) object@fit)
setGeneric("fit<-", function(x, value) standardGeneric("fit<-"))
setReplaceMethod("fit", "ksvm", function(x, value) {
  x@fit <- value
  x
})

if(!isGeneric("lev")){
  if (is.function("lev"))
    fun <- lev
  else fun <- function(object) standardGeneric("lev")
  setGeneric("lev", fun)
}
setMethod("lev", "ksvm", function(object) object@lev)
setGeneric("lev<-", function(x, value) standardGeneric("lev<-"))
setReplaceMethod("lev", "ksvm", function(x, value) {
  x@lev <- value
  x
})

if(!isGeneric("prior")){
  if (is.function("prior"))
    fun <- prior
  else fun <- function(object) standardGeneric("prior")
  setGeneric("prior", fun)
}
setMethod("prior", "ksvm", function(object) object@prior)
setGeneric("prior<-", function(x, value) standardGeneric("prior<-"))
setReplaceMethod("prior", "ksvm", function(x, value) {
  x@prior <- value
  x
})


if(!isGeneric("nclass")){
  if (is.function("nclass"))
    fun <- nclass
  else fun <- function(object) standardGeneric("nclass")
  setGeneric("nclass", fun)
}
setMethod("nclass", "ksvm", function(object) object@nclass)
setGeneric("nclass<-", function(x, value) standardGeneric("nclass<-"))
setReplaceMethod("nclass", "ksvm", function(x, value) {
  x@nclass <- value
  x
})

if(!isGeneric("alpha")){
  if (is.function("alpha"))
    fun <- alpha
  else fun <- function(object) standardGeneric("alpha")
  setGeneric("alpha", fun)
}
setMethod("alpha", "ksvm", function(object) object@alpha)
setGeneric("alpha<-", function(x, value) standardGeneric("alpha<-"))
setReplaceMethod("alpha", "ksvm", function(x, value) {
  x@alpha <- value
  x
})

if(!isGeneric("coeff")){
  if (is.function("coeff"))
    fun <- coeff
  else fun <- function(object) standardGeneric("coeff")
  setGeneric("coeff", fun)
}
setMethod("coeff", "ksvm", function(object) object@coeff)
setGeneric("coeff<-", function(x, value) standardGeneric("coeff<-"))
setReplaceMethod("coeff", "ksvm", function(x, value) {
  x@coeff <- value
  x
})

if(!isGeneric("alphaindex")){
  if (is.function("alphaindex"))
    fun <- alphaindex
  else fun <- function(object) standardGeneric("alphaindex")
  setGeneric("alphaindex", fun)
}
setMethod("alphaindex", "ksvm", function(object) object@alphaindex)
setGeneric("alphaindex<-", function(x, value) standardGeneric("alphaindex<-"))
setReplaceMethod("alphaindex", "ksvm", function(x, value) {
  x@alphaindex <- value
  x
})

if(!isGeneric("prob.model")){
  if (is.function("prob.model"))
    fun <- prob.model
  else fun <- function(object) standardGeneric("prob.model")
  setGeneric("prob.model", fun)
}
setMethod("prob.model", "ksvm", function(object) object@prob.model)
setGeneric("prob.model<-", function(x, value) standardGeneric("prob.model<-"))
setReplaceMethod("prob.model", "ksvm", function(x, value) {
  x@prob.model <- value
  x
})


if(!isGeneric("b")){
  if (is.function("b"))
    fun <- b
  else fun <- function(object) standardGeneric("b")
  setGeneric("b", fun)
}
setMethod("b", "ksvm", function(object) object@b)
setGeneric("b<-", function(x, value) standardGeneric("b<-"))
setReplaceMethod("b", "ksvm", function(x, value) {
  x@b <- value
  x
})

if(!isGeneric("SVindex")){
  if (is.function("SVindex"))
    fun <- SVindex
  else fun <- function(object) standardGeneric("SVindex")
  setGeneric("SVindex", fun)
}
setMethod("SVindex", "ksvm", function(object) object@SVindex)
setGeneric("SVindex<-", function(x, value) standardGeneric("SVindex<-"))
setReplaceMethod("SVindex", "ksvm", function(x, value) {
  x@SVindex <- value
  x
})

if(!isGeneric("nSV")){
  if (is.function("nSV"))
    fun <- nSV
  else fun <- function(object) standardGeneric("nSV")
  setGeneric("nSV", fun)
}
setMethod("nSV", "ksvm", function(object) object@nSV)
setGeneric("nSV<-", function(x, value) standardGeneric("nSV<-"))
setReplaceMethod("nSV", "ksvm", function(x, value) {
  x@nSV <- value
  x
})

if(!isGeneric("error")){
  if (is.function("error"))
    fun <- error
  else fun <- function(object) standardGeneric("error")
  setGeneric("error", fun)
}
setMethod("error", "ksvm", function(object) object@error)
setGeneric("error<-", function(x, value) standardGeneric("error<-"))
setReplaceMethod("error", "ksvm", function(x, value) {
  x@error <- value
  x
})

if(!isGeneric("cross")){
  if (is.function("cross"))
    fun <- cross
  else fun <- function(object) standardGeneric("cross")
  setGeneric("cross", fun)
}
setMethod("cross", "ksvm", function(object) object@cross)
setGeneric("cross<-", function(x, value) standardGeneric("cross<-"))
setReplaceMethod("cross", "ksvm", function(x, value) {
  x@cross <- value
  x
})

if(!isGeneric("n.action")){
  if (is.function("n.action"))
    fun <- n.action
  else fun <- function(object) standardGeneric("n.action")
  setGeneric("n.action", fun)
}
setMethod("n.action", "ksvm", function(object) object@n.action)
setGeneric("n.action<-", function(x, value) standardGeneric("n.action<-"))
setReplaceMethod("n.action", "ksvm", function(x, value) {
  x@n.action <- value
  x
})













## failed attempt to get rid of all this above


## mkaccesfun <- function(cls)
#{
#  snames <- slotNames(cls)
## 
#
#  for(i in 1:length(snames))
#    { resF <- paste("\"",snames[i],"\"",sep="")
#      if(!isGeneric(snames[i]))
#        eval(parse(file="",text=paste("setGeneric(",resF,",function(object)","standardGeneric(",resF,")",")",sep=" ")))   
#    setGeneric(snames[i], function(object) standardGeneric(snames[i]))
# 
#  setMethod(snames[i], cls, function(object)  eval(parse(file="",text=paste("object@",snames[i],sep=""))))
#  resG  <-  paste("\"",snames[i],"<-","\"",sep="")
#eval(parse(file="",text=paste("setGeneric(",resG,",function(x, value)","standardGeneric(",resG,")",")",sep=" ")))
#  setReplaceMethod(snames[i], cls, function(x, value) {
#    eval(parse(file="",text=paste("x@",snames[i],"<-value",sep="")))
#    x
#  })                   
#  }
#}


#kernel principal components object
setClass("kpca", representation(pcv = "matrix",
                                eig = "vector",
                                kernelf = "function",
                                kpar = "list",
                                rotated = "matrix",
                                xmatrix = "matrix",
                                kcall = "ANY",
                                kterms = "ANY",
                                n.action = "ANY"))
#accessor functions 
if(!isGeneric("pcv")){
  if (is.function("pcv"))
    fun <- pcv
  else fun <- function(object) standardGeneric("pcv")
  setGeneric("pcv", fun)
}
setMethod("pcv", "kpca", function(object) object@pcv)
setGeneric("pcv<-", function(x, value) standardGeneric("pcv<-"))
setReplaceMethod("pcv", "kpca", function(x, value) {
  x@pcv <- value
  x
})

if(!isGeneric("eig")){
  if (is.function("eig"))
    fun <- eig
  else fun <- function(object) standardGeneric("eig")
  setGeneric("eig", fun)
}
setMethod("eig", "kpca", function(object) object@eig)
setGeneric("eig<-", function(x, value) standardGeneric("eig<-"))
setReplaceMethod("eig", "kpca", function(x, value) {
  x@eig <- value
  x
})

if(!isGeneric("rotated")){
  if (is.function("rotated"))
    fun <- rotated
  else fun <- function(object) standardGeneric("rotated")
  setGeneric("rotated", fun)
}
setMethod("rotated", "kpca", function(object) object@rotated)
setGeneric("rotated<-", function(x, value) standardGeneric("rotated<-"))
setReplaceMethod("rotated", "kpca", function(x, value) {
  x@rotated <- value
  x
})



setMethod("kernelf","kpca", function(object) object@kernelf)
setReplaceMethod("kernelf","kpca", function(x, value){
  x@kernelf <- value
  x
})

setMethod("xmatrix","kpca", function(object) object@xmatrix)
setReplaceMethod("xmatrix","kpca", function(x, value){
  x@xmatrix <- value
  x
})

setMethod("kcall","kpca", function(object) object@kcall)
setReplaceMethod("kcall","kpca", function(x, value){
  x@kcall <- value
  x
})

setMethod("kterms","kpca", function(object) object@kterms)
setReplaceMethod("kterms","kpca", function(x, value){
  x@kterms <- value
  x
})

setMethod("n.action","kpca", function(object) object@n.action)
setReplaceMethod("n.action","kpca", function(x, value){
  x@n.action <- value
  x
})


setClass("ipop", representation(primal = "vector",
                                dual = "numeric",
                                how = "character"
                               ))

if(!isGeneric("primal")){
  if (is.function("primal"))
    fun <- primal
  else fun <- function(object) standardGeneric("primal")
  setGeneric("primal", fun)
}
setMethod("primal", "ipop", function(object) object@primal)
setGeneric("primal<-", function(x, value) standardGeneric("primal<-"))
setReplaceMethod("primal", "ipop", function(x, value) {
  x@primal <- value
  x
})

if(!isGeneric("dual")){
  if (is.function("dual"))
    fun <- dual
  else fun <- function(object) standardGeneric("dual")
  setGeneric("dual", fun)
}
setMethod("dual", "ipop", function(object) object@dual)
setGeneric("dual<-", function(x, value) standardGeneric("dual<-"))
setReplaceMethod("dual", "ipop", function(x, value) {
  x@dual <- value
  x
})

if(!isGeneric("how")){
  if (is.function("how"))
    fun <- how
  else fun <- function(object) standardGeneric("how")
  setGeneric("how", fun)
}
setMethod("how", "ipop", function(object) object@how)
setGeneric("how<-", function(x, value) standardGeneric("how<-"))
setReplaceMethod("how", "ipop", function(x, value) {
  x@how <- value
  x
})

# Kernel Canonical Correlation Analysis
setClass("kcca", representation(kcor = "vector",
                                xcoef = "matrix",
                                ycoef = "matrix",
                                xvar = "matrix",
                                yvar = "matrix"))


if(!isGeneric("kcor")){
  if (is.function("kcor"))
    fun <- kcor
  else fun <- function(object) standardGeneric("kcor")
  setGeneric("kcor", fun)
}
setMethod("kcor", "kcca", function(object) object@kcor)
setGeneric("kcor<-", function(x, value) standardGeneric("kcor<-"))
setReplaceMethod("kcor", "kcca", function(x, value) {
  x@kcor <- value
  x
})

if(!isGeneric("xcoef")){
  if (is.function("xcoef"))
    fun <- xcoef
  else fun <- function(object) standardGeneric("xcoef")
  setGeneric("xcoef", fun)
}
setMethod("xcoef", "kcca", function(object) object@xcoef)
setGeneric("xcoef<-", function(x, value) standardGeneric("xcoef<-"))
setReplaceMethod("xcoef", "kcca", function(x, value) {
  x@xcoef <- value
  x
})

if(!isGeneric("ycoef")){
  if (is.function("ycoef"))
    fun <- ycoef
  else fun <- function(object) standardGeneric("ycoef")
  setGeneric("ycoef", fun)
}
setMethod("ycoef", "kcca", function(object) object@ycoef)
setGeneric("ycoef<-", function(x, value) standardGeneric("ycoef<-"))
setReplaceMethod("ycoef", "kcca", function(x, value) {
  x@ycoef <- value
  x
})

if(!isGeneric("xvar")){
  if (is.function("xvar"))
    fun <- xvar
  else fun <- function(object) standardGeneric("xvar")
  setGeneric("xvar", fun)
}
setMethod("xvar", "kcca", function(object) object@xvar)
setGeneric("xvar<-", function(x, value) standardGeneric("xvar<-"))
setReplaceMethod("xvar", "kcca", function(x, value) {
  x@xvar <- value
  x
})

if(!isGeneric("yvar")){
  if (is.function("yvar"))
    fun <- yvar
  else fun <- function(object) standardGeneric("yvar")
  setGeneric("yvar", fun)
}
setMethod("yvar", "kcca", function(object) object@yvar)
setGeneric("yvar<-", function(x, value) standardGeneric("yvar<-"))
setReplaceMethod("yvar", "kcca", function(x, value) {
  x@yvar <- value
  x
})


setClass("gausspr",representation(tol = "numeric",
                                  kernelf = "function",
                                  kpar = "list",
                                  kcall = "ANY",
                                  type = "character",
                                  kterms = "ANY",
                                  xmatrix = "matrix",
                                  ymatrix = "ANY",
                                  fit = "ANY",
                                  lev = "vector",
                                  nclass = "numeric",
                                  alpha = "ANY",
                                  alphaindex="list",
                                  nvar = "numeric",
                                  error = "numeric",
                                  cross = "numeric",
                                  n.action = "ANY"))


         
setMethod("kpar","gausspr", function(object) object@kpar)
setReplaceMethod("kpar","gausspr", function(x, value){
  x@kpar <- value
  x
})

setMethod("lev","gausspr", function(object) object@lev)
setReplaceMethod("lev","gausspr", function(x, value){
  x@lev <- value
  x
})

setMethod("type","gausspr", function(object) object@type)
setReplaceMethod("type","gausspr", function(x, value){
  x@type <- value
  x
})

setMethod("kernelf","gausspr", function(object) object@kernelf)
setReplaceMethod("kernelf","gausspr", function(x, value){
  x@kernelf <- value
  x
})
         
setMethod("alpha","gausspr", function(object) object@alpha)
setReplaceMethod("alpha","gausspr", function(x, value){
  x@alpha <- value
  x
})

setMethod("xmatrix","gausspr", function(object) object@xmatrix)
setReplaceMethod("xmatrix","gausspr", function(x, value){
  x@xmatrix <- value
  x
})

setMethod("ymatrix","gausspr", function(object) object@ymatrix)
setReplaceMethod("ymatrix","gausspr", function(x, value){
  x@ymatrix <- value
  x
})

setMethod("nclass","gausspr", function(object) object@nclass)
setReplaceMethod("nclass","gausspr", function(x, value){
  x@nclass <- value
  x
})


setMethod("kcall","gausspr", function(object) object@kcall)
setReplaceMethod("kcall","gausspr", function(x, value){
  x@kcall <- value
  x
})

setMethod("fit","gausspr", function(object) object@fit)
setReplaceMethod("fit","gausspr", function(x, value){
  x@fit <- value
  x
})

setMethod("error","gausspr", function(object) object@error)
setReplaceMethod("error","gausspr", function(x, value){
  x@error <- value
  x
})

setMethod("cross","gausspr", function(object) object@cross)
setReplaceMethod("cross","gausspr", function(x, value){
  x@cross <- value
  x
})

setMethod("alphaindex","gausspr", function(object) object@alphaindex)
setReplaceMethod("alphaindex","gausspr", function(x, value){
  x@alphaindex <- value
  x
})



setMethod("kterms","gausspr", function(object) object@kterms)
setReplaceMethod("kterms","gausspr", function(x, value){
  x@kterms <- value
  x
})

setMethod("n.action","gausspr", function(object) object@n.action)
setReplaceMethod("n.action","gausspr", function(x, value){
  x@n.action <- value
  x
})         

# Relevance Vector Machine object 
setClass("rvm", representation( 
                                tol = "numeric",
                                kernelf = "function",
                                kpar = "list",
                                kcall = "ANY",
                                type = "character",
                                kterms = "ANY",
                                xmatrix = "matrix",
                                ymatrix = "ANY",
                                fit = "ANY",
                                lev = "vector",
                                nclass = "numeric",
                                alpha = "ANY",
                                nvar = "numeric",
                                mlike = "numeric",
                                RVindex = "vector",
                                nRV = "numeric",
                                cross = "ANY",
                                error = "numeric",
                                n.action= "ANY"))


if(!isGeneric("tol")){
  if (is.function("tol"))
    fun <- tol
  else fun <- function(object) standardGeneric("tol")
  setGeneric("tol", fun)
}
setMethod("tol", "rvm", function(object) object@tol)
setGeneric("tol<-", function(x, value) standardGeneric("tol<-"))
setReplaceMethod("tol", "rvm", function(x, value) {
  x@tol <- value
  x
})


setMethod("type","rvm", function(object) object@type)
setReplaceMethod("type","rvm", function(x, value){
  x@type <- value
  x
})

if(!isGeneric("RVindex")){
  if (is.function("RVindex"))
    fun <- RVindex
  else fun <- function(object) standardGeneric("RVindex")
  setGeneric("RVindex", fun)
}
setMethod("RVindex", "rvm", function(object) object@RVindex)
setGeneric("RVindex<-", function(x, value) standardGeneric("RVindex<-"))
setReplaceMethod("RVindex", "rvm", function(x, value) {
  x@RVindex <- value
  x
})

if(!isGeneric("nvar")){
  if (is.function("nvar"))
    fun <- nvar
  else fun <- function(object) standardGeneric("nvar")
  setGeneric("nvar", fun)
}
setMethod("nvar", "rvm", function(object) object@nvar)
setGeneric("nvar<-", function(x, value) standardGeneric("nvar<-"))
setReplaceMethod("nvar", "rvm", function(x, value) {
  x@nvar <- value
  x
})

if(!isGeneric("nRV")){
  if (is.function("nRV"))
    fun <- nRV
  else fun <- function(object) standardGeneric("nRV")
  setGeneric("nRV", fun)
}
setMethod("nRV", "rvm", function(object) object@nRV)
setGeneric("nRV<-", function(x, value) standardGeneric("nRV<-"))
setReplaceMethod("nRV", "rvm", function(x, value) {
  x@nRV <- value
  x
})


if(!isGeneric("mlike")){
  if (is.function("mlike"))
    fun <- mlike
  else fun <- function(object) standardGeneric("mlike")
  setGeneric("mlike", fun)
}
setMethod("mlike", "rvm", function(object) object@mlike)
setGeneric("mlike<-", function(x, value) standardGeneric("mlike<-"))
setReplaceMethod("mlike", "rvm", function(x, value) {
  x@mlike <- value
  x
})


setMethod("kpar","rvm", function(object) object@kpar)
setReplaceMethod("kpar","rvm", function(x, value){
  x@kpar <- value
  x
})

setMethod("lev","rvm", function(object) object@lev)
setReplaceMethod("lev","rvm", function(x, value){
  x@lev <- value
  x
})

setMethod("kernelf","rvm", function(object) object@kernelf)
setReplaceMethod("kernelf","rvm", function(x, value){
  x@kernelf <- value
  x
})

setMethod("xmatrix","rvm", function(object) object@xmatrix)
setReplaceMethod("xmatrix","rvm", function(x, value){
  x@xmatrix <- value
  x
})

setMethod("ymatrix","rvm", function(object) object@ymatrix)
setReplaceMethod("ymatrix","rvm", function(x, value){
  x@ymatrix <- value
  x
})

setMethod("nclass","rvm", function(object) object@nclass)
setReplaceMethod("nclass","rvm", function(x, value){
  x@nclass <- value
  x
})

setMethod("kcall","rvm", function(object) object@kcall)
setReplaceMethod("kcall","rvm", function(x, value){
  x@kcall <- value
  x
})

setMethod("fit","rvm", function(object) object@fit)
setReplaceMethod("fit","rvm", function(x, value){
  x@fit <- value
  x
})

setMethod("kterms","rvm", function(object) object@kterms)
setReplaceMethod("kterms","rvm", function(x, value){
  x@kterms <- value
  x
})
                 

setMethod("alpha","rvm", function(object) object@alpha)
setReplaceMethod("alpha","rvm", function(x, value){
  x@alpha <- value
  x
})


setMethod("error","rvm", function(object) object@error)
setReplaceMethod("error","rvm", function(x, value){
  x@error <- value
  x
})

setMethod("cross","rvm", function(object) object@cross)
setReplaceMethod("cross","rvm", function(x, value){
  x@cross <- value
  x
})

setMethod("n.action","rvm", function(object) object@n.action)
setReplaceMethod("n.action","rvm", function(x, value){
  x@n.action <- value
  x
})

setClass("inc.chol",representation("matrix",
 				    pivots="vector",
				    diag.residues="vector",
				    maxresiduals="vector"),
				    prototype=structure(.Data=matrix(),
				    pivots=vector(),
				    diag.residues=vector(), 
				    maxresiduals=vector()))


if(!isGeneric("pivots")){
if (is.function("pivots"))
  fun <- pivots
else fun <- function(object) standardGeneric("pivots")
setGeneric("pivots", fun)
}
setMethod("pivots", "inc.chol", function(object) object@pivots)
setGeneric("pivots<-", function(x, value) standardGeneric("pivots<-"))
setReplaceMethod("pivots", "inc.chol", function(x, value) {
  x@pivots <- value
  x
})

if(!isGeneric("diag.residues")){
if (is.function("diag.residues"))
  fun <- diag.residues
else fun <- function(object) standardGeneric("diag.residues")
setGeneric("diag.residues", fun)
}
setMethod("diag.residues", "inc.chol", function(object) object@diag.residues)
setGeneric("diag.residues<-", function(x,value) standardGeneric("diag.residues<-"))
setReplaceMethod("diag.residues", "inc.chol", function(x, value) {
  x@diag.residues <- value
  x
})

if(!isGeneric("maxresiduals")){
if (is.function("maxresiduals"))
  fun <- maxresiduals
else fun <- function(object) standardGeneric("maxresiduals")
setGeneric("maxresiduals", fun)
}
setMethod("maxresiduals", "inc.chol", function(object) object@maxresiduals)
setGeneric("maxresiduals<-", function(x,value) standardGeneric("maxresiduals<-"))
setReplaceMethod("maxresiduals", "inc.chol", function(x, value) {
  x@maxresiduals <- value
  x
})




setClass("specc",representation("vector",
                                centers="matrix",
                                size="vector",
                                kernelf="function",
                                withinss = "vector"
                                ),prototype=structure(.Data=vector(),
                                    centers = matrix(),
                                    size=matrix(),
                                    kernelf = ls,
                                    withinss=vector()))


if(!isGeneric("centers")){
if (is.function("centers"))
  fun <- centers
else fun <- function(object) standardGeneric("centers")
setGeneric("centers", fun)
}
setMethod("centers", "specc", function(object) object@centers)
setGeneric("centers<-", function(x,value) standardGeneric("centers<-"))
setReplaceMethod("centers", "specc", function(x, value) {
  x@centers <- value
  x
})

if(!isGeneric("size")){
if (is.function("size"))
  fun <- size
else fun <- function(object) standardGeneric("size")
setGeneric("size", fun)
}
setMethod("size", "specc", function(object) object@size)
setGeneric("size<-", function(x,value) standardGeneric("size<-"))
setReplaceMethod("size", "specc", function(x, value) {
  x@size <- value
  x
})

if(!isGeneric("withinss")){
if (is.function("withinss"))
  fun <- withinss
else fun <- function(object) standardGeneric("withinss")
setGeneric("withinss", fun)
}
setMethod("withinss", "specc", function(object) object@withinss)
setGeneric("withinss<-", function(x,value) standardGeneric("withinss<-"))
setReplaceMethod("withinss", "specc", function(x, value) {
  x@withinss <- value
  x
})

setMethod("kernelf","specc", function(object) object@kernelf)
setReplaceMethod("kernelf","specc", function(x, value){
  x@kernelf <- value
  x
})



setClass("ranking",representation("matrix",
 				    convergence="matrix",
				    edgegraph="matrix"),
				    prototype=structure(.Data=matrix(),
                                    convergence=matrix(),
				    edgegraph=matrix()))

if(!isGeneric("convergence")){
if (is.function("convergence"))
  fun <- convergence
else fun <- function(object) standardGeneric("convergence")
setGeneric("convergence", fun)
}
setMethod("convergence", "ranking", function(object) object@convergence)
setGeneric("convergence<-", function(x,value) standardGeneric("convergence<-"))
setReplaceMethod("convergence", "ranking", function(x, value) {
  x@convergence <- value
  x
})

if(!isGeneric("edgegraph")){
if (is.function("edgegraph"))
  fun <- edgegraph
else fun <- function(object) standardGeneric("edgegraph")
setGeneric("edgegraph", fun)
}
setMethod("edgegraph", "ranking", function(object) object@edgegraph)
setGeneric("edgegraph<-", function(x,value) standardGeneric("edgegraph<-"))
setReplaceMethod("edgegraph", "ranking", function(x, value) {
  x@edgegraph <- value
  x
})

## online learning algorithms class

setClass("onlearn", representation(
                                kernelf = "function",
                                buffer = "numeric",
                                kpar = "list",
                                xmatrix = "matrix",
                                fit = "numeric",
                                onstart = "numeric",
                                onstop = "numeric",
                                alpha = "ANY",
                                rho = "numeric",
                                b = "numeric",
                                pattern ="factor",
                                type="character"
                               ))


setMethod("fit","onlearn", function(object) object@fit)
setReplaceMethod("fit","onlearn", function(x, value){
  x@fit <- value
  x
})

if(!isGeneric("onstart")){
  if (is.function("onstart"))
    fun <- onstart
  else fun <- function(object) standardGeneric("onstart")
  setGeneric("onstart", fun)
}
setMethod("onstart", "onlearn", function(object) object@onstart)
setGeneric("onstart<-", function(x, value) standardGeneric("onstart<-"))
setReplaceMethod("onstart", "onlearn", function(x, value) {
  x@onstart <- value
  x
})

if(!isGeneric("onstop")){
  if (is.function("onstop"))
    fun <- onstop
  else fun <- function(object) standardGeneric("onstop")
  setGeneric("onstop", fun)
}
setMethod("onstop", "onlearn", function(object) object@onstop)
setGeneric("onstop<-", function(x, value) standardGeneric("onstop<-"))
setReplaceMethod("onstop", "onlearn", function(x, value) {
  x@onstop <- value
  x
})

if(!isGeneric("buffer")){
  if (is.function("buffer"))
    fun <- buffer
  else fun <- function(object) standardGeneric("buffer")
  setGeneric("buffer", fun)
}
setMethod("buffer", "onlearn", function(object) object@buffer)
setGeneric("buffer<-", function(x, value) standardGeneric("buffer<-"))
setReplaceMethod("buffer", "onlearn", function(x, value) {
  x@buffer <- value
  x
})

setMethod("kernelf","onlearn", function(object) object@kernelf)
setReplaceMethod("kernelf","onlearn", function(x, value){
  x@kernelf <- value
  x
})

setMethod("kpar","onlearn", function(object) object@kpar)
setReplaceMethod("kpar","onlearn", function(x, value){
  x@kpar <- value
  x
})

setMethod("xmatrix","onlearn", function(object) object@xmatrix)
setReplaceMethod("xmatrix","onlearn", function(x, value){
  x@xmatrix <- value
  x
})


setMethod("alpha","onlearn", function(object) object@alpha)
setReplaceMethod("alpha","onlearn", function(x, value){
  x@alpha <- value
  x
})

setMethod("b","onlearn", function(object) object@b)
setReplaceMethod("b","onlearn", function(x, value){
  x@b <- value
  x
})

setMethod("type","onlearn", function(object) object@type)
setReplaceMethod("type","onlearn", function(x, value){
  x@type <- value
  x
})

if(!isGeneric("rho")){
  if (is.function("rho"))
    fun <- rho
  else fun <- function(object) standardGeneric("rho")
  setGeneric("rho", fun)
}
setMethod("rho", "onlearn", function(object) object@rho)
setGeneric("rho<-", function(x, value) standardGeneric("rho<-"))
setReplaceMethod("rho", "onlearn", function(x, value) {
  x@rho <- value
  x
})

if(!isGeneric("pattern")){
  if (is.function("pattern"))
    fun <- pattern
  else fun <- function(object) standardGeneric("pattern")
  setGeneric("pattern", fun)
}
setMethod("pattern", "onlearn", function(object) object@pattern)
setGeneric("pattern<-", function(x, value) standardGeneric("pattern<-"))
setReplaceMethod("pattern", "onlearn", function(x, value) {
  x@pattern <- value
  x
})





setClass("kfa",representation(alpha = "matrix",
                              alphaindex = "vector",
                              kernelf = "function",
                              xmatrix = "matrix",
                              kcall = "ANY",
                              kterms = "ANY" )) 


setMethod("kernelf","kfa", function(object) object@kernelf)
setReplaceMethod("kernelf","kfa", function(x, value){
  x@kernelf <- value
  x
})

setMethod("alphaindex","kfa", function(object) object@alphaindex)
setReplaceMethod("alphaindex","kfa", function(x, value){
  x@alphaindex <- value
  x
})

setMethod("alpha","kfa", function(object) object@alpha)
setReplaceMethod("alpha","kfa", function(x, value){
  x@alpha <- value
  x
})

setMethod("xmatrix","kfa", function(object) object@xmatrix)
setReplaceMethod("xmatrix","kfa", function(x, value){
  x@xmatrix <- value
  x
})


setMethod("kcall","kfa", function(object) object@kcall)
setReplaceMethod("kcall","kfa", function(x, value){
  x@kcall <- value
  x
})


setMethod("kterms","kfa", function(object) object@kterms)
setReplaceMethod("kterms","kfa", function(x, value){
  x@kterms <- value
  x
})
