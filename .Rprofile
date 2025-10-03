# Error/Warning colorization for visibility
options(
  error = function() {
    msg <- geterrmessage()
    cat("\033[31m", "ERROR: ", msg, "\033[0m\n", sep = "")
    invokeRestart("abort")
  }
)

options(
  warning.expression = quote({
    msg <- geterrmessage()
    cat("\033[33m", "WARNING: ", msg, "\033[0m\n", sep = "")
  })
)
