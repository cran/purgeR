# Unload dynamic libraries on unload
.onUnload <- function(libpath) {
  library.dynam.unload("purgeR", libpath)
}
