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

To run them in GitHub Actions on a pull request, either:

* add the `run-extended` label to the PR (re-triggers R-CMD-check with the 
  env var set), or
* manually dispatch the workflow with the `extended` input set to `true`, e.g.
  `gh workflow run R-CMD-check.yaml --ref <branch> -f extended=true`.
