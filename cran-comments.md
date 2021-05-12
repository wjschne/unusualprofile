## Test environments

-   local OS X install, R 3.5.1
-   ubuntu 12.04 (on travis-ci), R 3.5.1
-   win-builder (devel and release)

## R CMD check results

Duration: 36.4s

0 errors v \| 0 warnings v \| 0 notes v

R CMD check succeeded

## Addressing concerns from a previous build

From Gregor Seyer:

> Please add () behind all function names in the description texts 
(DESCRIPTION file). e.g: --> cond_maha()

Done

> Please omit the space within the doi specification to make it clickable.

Done

> Please add \value to .Rd files regarding exported methods and explain 
the functions results in the documentation. Please write about the 
structure of the output (class) and also what the output means. (If a 
function does not return a value, please document that too, e.g. 
\value{No return value, called for side effects} or similar)
Missing Rd-tags:
      plot.cond_maha.Rd: \value
      plot.maha.Rd: \value


Rd-tags were completed for both functions

> Please always make sure to reset to user's options(), working directory 
or par() after you changed it in examples and vignettes and demos.
e.g.: inst/doc/..
old <- options(digits = 3)
...
options(old)

I removed the options() command.

