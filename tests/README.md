### Running extended tests

The testthat directory contains several long-running tests which are not run by 
default. Running all of these tests can take several hours on an ordinary 
laptop. These extended tests are gated on the environment variable 
`COEVOLVE_EXTENDED_TESTS = "true"`.

To run them locally:

```r
Sys.setenv(COEVOLVE_EXTENDED_TESTS = "true")
devtools::test()
```

In GitHub Actions, the extended tests run as a separate `extended-tests`
workflow that splits each long-running test file into its own parallel job.
Trigger it by either:

* adding the `run-extended` label to the PR, or
* dispatching the workflow manually:
  `gh workflow run extended-tests.yaml --ref <branch>`.

R-CMD-check itself never runs the extended tests, so the standard CI signal
stays fast regardless of the label.
