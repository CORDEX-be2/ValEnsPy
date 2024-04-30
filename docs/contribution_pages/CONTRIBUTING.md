# Contributing

All contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

## Ways to contribute
There are many ways to contribute to the ValEnsPy package:
1. Partake in one of the [discussions](https://github.com/CORDEX-be2/ValEnsPy/discussions) or create a new one.
2. Create a new [issue](https://github.com/CORDEX-be2/ValEnsPy/issues) to help us improve the package.
  - New feature requests
  - Bug Reports

3. Contribute to the codebase. This can be done by:
  - Taking on an open issue
  - Adding a new diagnostic (see ???)

4. Help improve the documentation. This can be done by:
- Adding a new example
- Improving the general documentation

## Contribute to the codebase
Ready to make code contributions? Here's how to set up ValEnsPy for local development.
If you are planning to develop on the VSC see CONTRIBUTING_VSC.md for a VSC specific guide.

### Required software

The following software (or equivalent) is required to set up a developer environment:
* [Anaconda](https://anaconda.org/)

Make sure you have this software installed before proceeding.

### Setup a developer environment
1. Clone the ValEnsPy locally:

```bash
git clone git@github.com:CORDEX-be2/ValEnsPy.git
```
2. Create a conda environment and install the required packages.
```bash
# Setup a developers' environment
conda create -n valenspy_dev python==3.9 poetry
conda activate valenspy_dev

# Install dependencies in the developers' environment
cd ValEnsPy
poetry install
```
3. Create a branch for local development which is a copy of the **dev** branch.:
```bash
# checkout the dev branch
git checkout dev
git pull

# Create a new local branch and switch to it.
git branch name-of-your-bugfix-or-feature
git checkout name-of-your-bugfix-or-feature
```
Now you can make local changes.

5. Push your code online:
```bash
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
From time to time the *dev* branch will be merged with the master with a new [*Release tag*](https://github.com/CORDEX-be2/ValEnsPy/releases).

## Acknowledgement
This file is inspired by the [RavenPy](https://github.com/CSHS-CWRA/RavenPy) project. Thank you for the inspiration!”.