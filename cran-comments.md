## Resubmission

This is a resubmission. In this version I have:

- Changed the Description field of the DESCRIPTION file to try and better explain the method of fitting
    - The package does not use a new method of non-linear least square regression. Instead it uses nlsLM which uses the Levenberg-Marquadt algorithm, so we do not think a reference is needed
    - We did however, expand on the Description field: "Non-linear least squares regression with the Levenberg-Maquardt algorithm using multiple starting values for increasing the chance that the minimum found is the global minimum."

## Test environments

- local OS X install, R 3.4.3
- unbuntu 14.04.5 (on travis-ci), R 3.4.2
- Mac OS X (on tavis-ci), R 3.4.3
- win-builder (using devtools::build_win())

## R CMD check results

- There were no ERRORs or WARNINGs
- There was 1 NOTE:
    - This is purely because it is a New Submission
