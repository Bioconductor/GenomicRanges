### =========================================================================
### Some low-level (non exported) utility functions.
### -------------------------------------------------------------------------


hasHead <- function(x, h) {
  identical(head(x, length(h)), h)
}

