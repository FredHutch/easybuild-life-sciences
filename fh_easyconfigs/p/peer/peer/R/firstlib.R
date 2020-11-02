.First.lib <- function(lib, pkg) {
  library.dynam("libpeer", pkg, lib)
}
