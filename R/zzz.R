## Startup code

.onUnload <- function(libpath) {
  library.dynam.unload("clustTMB", libpath)
}
