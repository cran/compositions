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
cat("Please be advised that this peace of software does complicated\n")
cat("geometric operations in high dimensional spaces often based on\n")
cat("on various incompatible geometries in most cases behaving more\n")
cat("like an expert system then like an application software merily\n")
cat("guessing the your intention using several incompatible systems\n")
cat("of thought. Thus please never trust the software. Always check.\n")

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
cat("Please be advised that this peace of software does complicated\n")
cat("geometric operations in high dimensional spaces often based on\n")
cat("on various incompatible geometries in most cases behaving more\n")
cat("like an expert system then like an application software merily\n")
cat("guessing the your intention using several incompatible systems\n")
cat("of thought. Thus please never trust the software. Always check.\n")
}

