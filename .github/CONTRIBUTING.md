# Contributing to coevolve

Development is a community effort, and we welcome participation.

## Code of Conduct

Please note that the coevolve package is released with a 
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). 

## Discussions

At <https://github.com/ScottClaessens/coevolve/discussions>, you can post 
general questions, brainstorm ideas, and ask for help with the coevolve package.

## Issues

<https://github.com/ropensci/targets/issues> is for maintenance tasks and 
feature requests.

## Development

External code contributions are extremely helpful in the right circumstances. 
Here are the recommended steps.

1. Prior to contribution, please propose your idea in a discussion thread so you
and the maintainer can define the intent and scope of your work.
2. [Fork the repository](https://help.github.com/articles/fork-a-repo/).
3. Follow the [GitHub flow](https://guides.github.com/introduction/flow/index.html)
to create a new branch, add commits, and open a pull request.
4. Discuss your code with the maintainer in the pull request thread.
5. If everything looks good, the maintainer will merge your code into the 
project.

Please also follow these additional guidelines.

* Respect the architecture and reasoning of the package.
* If possible, keep contributions small enough to easily review manually. It is
okay to split up your work into multiple pull requests.
* Format your code according to the 
[tidyverse style guide](https://style.tidyverse.org/) and check your formatting 
with the `lint_package()` function from the 
[`lintr`](https://github.com/jimhester/lintr) package.
* For new features or functionality, add tests in `tests`.
* Check code coverage with `covr::package_coverage()`. Automated tests should 
cover all the new or changed functionality in your pull request.
* Run overall package checks with `devtools::check()` and `goodpractice::gp()`
* Describe your contribution in the project's 
[`NEWS.md`](https://github.com/ropensci/targets/blob/main/NEWS.md) file. Be sure
to mention relevent GitHub issue numbers and your GitHub name.
* If you feel contribution is substantial enough for official author or 
contributor status, please add yourself to the `Authors@R` field of the 
[`DESCRIPTION`](https://github.com/ropensci/targets/blob/main/DESCRIPTION) file.
