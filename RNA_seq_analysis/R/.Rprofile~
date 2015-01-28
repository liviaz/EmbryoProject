## Example of .Rprofile

options(width=65, digits=5)
options(show.signif.stars=FALSE)
setHook(packageEvent("grDevices", "onLoad"),
        function(...) grDevices::ps.options(horizontal=FALSE))
set.seed(1234)

options(defaultPackages = c("utils","grDevices","graphics","stats","cummeRbund","DESeq2")) 

.First <- function() cat("
   Welcome to R!
 
")

.Last <- function()  cat("
   Goodbye!
 
")
