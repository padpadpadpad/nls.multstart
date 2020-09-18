## Resubmission

This is a minor resubmission. In this version we have:

- updated the argument supp_errors so that warning messages are suppressed as well as error messages

From Uwe Ligges comments:
- Please change http --> https, add trailing slashes, or follow moved
content as appropriate.
    - changed http --> https in README.md
- Is there some reference about the method you can add in the Description field in the form Authors (year) <doi:.....>?
    - Not really. This method is being used in primary research but no paper exists to cite its use more generally. The R package currently has 17 citations and is in general use. Consequently there is no obvious reference to add for the method implemented.

## Test environments

- local OS X install, R 4.0.2
- ubuntu 16.04.6 (on travis-ci), R 4.0.2
- win-builder (using devtools::check_win_devel())

## R CMD check resultsjust 

- There were no ERRORs or WARNINGs
