## First resubmission (May 12, 2021)
We have addressed the two concerns mentioned in the mail we got on May 11, 2021 (thanks Gregor Seyer). Specifically,
* A more elaborate package description was added in the description field of the DESCRIPTION file. Some references were included.
* We have unwraped (removed \dontrun) the examples that take < 5 seconds to be executed. We have replaced \dontrun with \dontest for an example that takes > 5 seconds to be executed and included a short test (using \dontshow) so this function (modelcompK) can be automatically checked.

## Test envirnoments
* local Windows 10 x64, R 4.0.5
* Ubuntu Linux 16.04 LTS, R-release, GCC (check_rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (check_rhub)
* Fedora Linux, R-devel, clang, gfortran (check_rhub)
* (check_win_devel)

## R CMD check results
* There were no ERRORs or WARNINGs.
* There was 1 NOTE:
```
Mantainer: 'Pablo Nájera pablo.najera@uam.es'

New submission

Possibly mis-spelled words in DESCRIPTION: 
CDM (7:70, 7:265)
```
All possibly mis-spelled words are surnames, acronyms, or British spelling used in paper references.
