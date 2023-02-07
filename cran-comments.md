## First Submission to CRAN, version 0.0.1

The first submission of the `CoxAIPW` package

### Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

### R CMD check results
> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
* checking CRAN incoming feasibility ... NOTE

Maintainer: 'Jiyu Luo <charlesluo1002@gmail.com>'
New submission
 
 > On fedora-clang-devel (r-devel)
 * checking for detritus in the temp directory ... NOTE
 
Found the following files/directories:
  'lastMiKTeXException'

0 errors ✓ | 0 warnings ✓ | 2 notes x

Notes:
* The 'lastMiKTeXException' note seems to be due to a MiKTeX bug during checking and can probably be ignored. See similar issues: https://github.com/r-hub/rhub/issues/503 .
* In the email from the Ubuntu test it actually says PREPERROR, but once you click in the reports, the status is success. This is likely due to the logs being too long and are automatically truncated, so it probably can be ignored. See similar issues: https://github.com/r-hub/rhub/issues/448 .