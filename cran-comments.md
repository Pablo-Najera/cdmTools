## cdmTools_1.0.3 (December 15, 2022)

#### News & Comments
* We have slightly modified the CA.MI function
* We have added the RDINA function
* We have added internal functions

#### Test envirnoments
* Local Windows 10 x64, R 4.1.2
* Ubuntu Linux 16.04 LTS, R-release, GCC (check_rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (check_rhub)
* Fedora Linux, R-devel, clang, gfortran (check_rhub)
* (check_win_devel)

#### R CMD check results

* There were no ERRORs or WARNINGs.
* There were 4 NOTEs:
```
Maintainer: 'Pablo Nájera <pablo.najera@uam.es>'
   
Possibly misspelled words in DESCRIPTION:
Abad (7:213, 7:497)
CDM (7:70, 7:412, 7:595, 7:904)
Chiu (7:670)
Kreitchmann (7:803)
Nájera (7:182, 7:489)
Torre (7:204)
al (7:678, 7:818)
de (7:198)
et (7:675, 7:815)

Uses the superseded package: ‘doSNOW’
   
Found the following (possibly) invalid DOIs:
  DOI: 10.1111/bmsp.12228
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503
```
All possibly misspelled words are surnames, acronyms, or British spelling used in paper references. The DOI 10.1111/bmsp.12228 is valid (from the British Journal of Mathematical and Statistical Psychology). The `doSNOW` package is required to show a progress bar inside a foreach loop.

## cdmTools_1.0.2 (May 17, 2022)

#### News & Comments
* We have added the CA.MI function and Rodrigo S. Kreitchmann as a contributor
* {CDM} package dependeny has been removed - The {cdmTools} package was removed from CRAN due to an issue with the {CDM} package

#### Test envirnoments
* Local Windows 10 x64, R 4.1.2
* Ubuntu Linux 16.04 LTS, R-release, GCC (check_rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (check_rhub)
* Fedora Linux, R-devel, clang, gfortran (check_rhub)
* (check_win_devel)

#### R CMD check results

* There were no ERRORs or WARNINGs.
* There was 1 NOTE:
```
Maintainer: 'Pablo Nájera <pablo.najera@uam.es>'
   
Possibly misspelled words in DESCRIPTION:
Abad (7:213, 7:497)
CDM (7:70, 7:412, 7:595, 7:904)
Chiu (7:670)
Kreitchmann (7:803)
Nájera (7:182, 7:489)
Torre (7:204)
al (7:678, 7:818)
de (7:198)
et (7:675, 7:815)

Uses the superseded package: ‘doSNOW’
   
Found the following (possibly) invalid DOIs:
  DOI: 10.1111/bmsp.12228
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503
```
All possibly misspelled words are surnames, acronyms, or British spelling used in paper references. The DOI 10.1111/bmsp.12228 is valid (from the British Journal of Mathematical and Statistical Psychology). The `doSNOW` package is required to show a progress bar inside a foreach loop.

## cdmTools_1.0.1 (March 23, 2022)

#### NEWS
* We have added the GNPC function and print methods for classes used in the package

#### Test envirnoments
* Local Windows 10 x64, R 4.1.2
* Ubuntu Linux 16.04 LTS, R-release, GCC (check_rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (check_rhub)
* Fedora Linux, R-devel, clang, gfortran (check_rhub)
* (check_win_devel)

#### R CMD check results

* There were no ERRORs or WARNINGs.
* There was 1 NOTE:
```
Maintainer: 'Pablo Nájera <pablo.najera@uam.es>'
   
   Found the following (possibly) invalid DOIs:
     DOI: 10.1111/bmsp.12228
       From: DESCRIPTION
       Status: Service Unavailable
       Message: 503
```
The DOI 10.1111/bmsp.12228 is valid (from the British Journal of Mathematical and Statistical Psychology).

## cdmTools_1.0.0 (May 12, 2021)

#### First resubmission (May 12, 2021)
We have addressed the two concerns mentioned in the mail we got on May 11, 2021 (thanks Gregor Seyer). Specifically,
* A more elaborate package description was added in the description field of the DESCRIPTION file. Some references were included.
* We have unwraped (removed \dontrun) the examples that take < 5 seconds to be executed. We have replaced \dontrun with \dontest for an example that takes > 5 seconds to be executed and included a short test (using \dontshow) so this function (modelcompK) can be automatically checked.

#### Test envirnoments
* Local Windows 10 x64, R 4.0.5
* Ubuntu Linux 16.04 LTS, R-release, GCC (check_rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (check_rhub)
* Fedora Linux, R-devel, clang, gfortran (check_rhub)
* (check_win_devel)

#### R CMD check results
* There were no ERRORs or WARNINGs.
* There was 1 NOTE:
```
Mantainer: 'Pablo Nájera pablo.najera@uam.es'

New submission

Possibly mis-spelled words in DESCRIPTION: 
CDM (7:70, 7:265)
```
All possibly mis-spelled words are surnames, acronyms, or British spelling used in paper references.
