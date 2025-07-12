### Running extended tests

The testthat directory contains several long-running tests which are not run by 
default. Running all of these tests can take several hours on an ordinary 
laptop. These extended tests can be switched on by defining an environmental 
variable `COEVOLVE_EXTENDED_TESTS = "true"` (e.g., 
`Sys.setenv("COEVOLVE_EXTENDED_TESTS" = "true")` in R), and on GitHub Actions by 
adding `run-extended` to the commit message.
