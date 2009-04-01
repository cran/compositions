.First.lib <- function(libname,pkgname,...) {
if( is.null(getOption("robust"))) options(robust=FALSE)
cat("Welcome to compositions, a package for compositional data analysis.\n")
cat("\n");
cat("Get help with \"? compositions\"")
cat("\n");
cat("Version info: to ensure compatibility with R > 2.7.0.\n")
cat(" \".inv\" in function names has globaly been replaced by \"Inv\"\n")
cat(" similarly \".row\" and \".col\" respectively by \"Row\" and \"Col\".\n")
cat("\n")
cat("WARNING\n")
cat("Please, take into account that this piece of software does complex   \n")
cat("geometric operations in high-dimensional spaces, possibly using      \n")
cat("different incompatible geometries. In particular, be advised that    \n")
cat("arithmetic operations (+ and *) might have non-conventional meanings.\n")
cat("Double-check that you are doing reasonable things. The mere fact that\n")
cat("the package computes something does not imply that this is reasonable\n")

}

.onAttach <- function(libname,pkgname,...) {
if( is.null(getOption("robust"))) options(robust=FALSE)
lapply(objects(pattern="^gsi.*",envir=parent.env(environment())),function(x) if( ! (x %in% objects(pattern="^gsi","package:compositions"))) assign(x,get(x),gsi))
cat("Welcome to compositions, a package for compositional data analysis.\n")
cat("\n");
cat("Get help with \"? compositions\"")
cat("\n");
cat("Version info: to ensure compatibility with R > 2.7.0.\n")
cat(" \".inv\" in function names has globaly been replaced by \"Inv\"\n")
cat(" similarly \".row\" and \".col\" respectively by \"Row\" and \"Col\".\n")
cat("\n")
cat("WARNING\n")
cat("Please, take into account that this piece of software does complex   \n")
cat("geometric operations in high-dimensional spaces, possibly using      \n")
cat("different incompatible geometries. In particular, be advised that    \n")
cat("arithmetic operations (+ and *) might have non-conventional meanings.\n")
cat("Double-check that you are doing reasonable things. The mere fact that\n")
cat("the package computes something does not imply that this is reasonable\n")
}

