# Contributing

All contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

## Get Started
Ready to make code contributions? Here is how to set up a developer's environment for the toolkit.

### Required software

The following software (or equivalent) is required to set up a developer environment:
* [Anaconda](https://anaconda.org/)

Make sure you have this software installed before proceeding.

### Setup a developer environment
1. Clone the MetObs-toolkit locally:

  ```
  git clone git@github.com:CORDEX-be2/ValEnsPy.git
  ```
2. Create a conda environment and install the required packages.
  ```
  # Setup a developers' environment
  conda create -n metobs_dev python==3.9 poetry
  conda activate metobs_dev

  # Install dependencies in the developers' environment
  cd MetObs_toolkit
  poetry install
  ```
3. Create a branch for local development which is a copy of the **dev** branch.:
 ```
  # checkout the dev branch
  git checkout dev
  git pull

  # Create a new local branch and switch to it.
  git branch name-of-your-bugfix-or-feature
  git checkout name-of-your-bugfix-or-feature
  ```
 Now you can make local changes.

5. Push your code online:
   ```
   # Add your changes to your commit
   git add -A
   # Write commit text
   git commit -m "Some text describing your code changes in this commit"
   # Push your branch online
   #only the first time:
   git push --set-upstream origin name-of-your-bugfix-or-feature
   #all other times
   git push
   ```

## Pull Request Guidelines
Once your branch has been *pushed* to github, you can create a *Pull request* in github. Make sure that you have **referred the corresponding issues** to the *Pull request*.
If your code adaptations are still *work-in-progress* add the ![Static Badge](https://img.shields.io/badge/WIP%20-%20%23A21079) label to it. For each push, github will perform a list of checks (package building, version control, functionality test, os-tests, documentation build test), in order to merge your contributions these tests must all be successful.

If your code is ready for review, you can add the ![Static Badge](https://img.shields.io/badge/Ready_for_Review%20-%20%230315E4) label to it.

After the code review, and all review marks are resolved, your contributions will be merged to the *dev* branch.

 ## Versioning/Tagging
 From time to time the *dev* branch will be merged with the master with a new [*Release tag*](https://github.com/vergauwenthomas/MetObs_toolkit/releases). The new release will be deployed to [PyPi index](https://pypi.org/project/MetObs-toolkit/) with the adequate versioning specified.

## Support
For general support or questions, you can refer them to @vergauwenthomas, or by mail to (thomas.vergauwen@meteo.be).

## Acknowledgement
This file is inspired by the [RavenPy](https://github.com/CSHS-CWRA/RavenPy) project. Thank you for the inspiration!”.