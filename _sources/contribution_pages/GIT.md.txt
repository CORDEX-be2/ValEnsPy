# Mini Git Guide
The purpose of this page is to provide links to resource for those who are new to git.

## Why git and GitHub?
Git is a version control system that allows you to track changes in your code. This allows you to work on a project with multiple people simultaneously without overwriting each other's work.

GitHub is a platform that allows you to store your code online. It provides a way to collaborate with others through issues, pull requests, and much more.

## Resources
An introduction to git
- [Intro to GIT](https://learngitbranching.js.org/)
A beginner friendly interactive tutorial that teaches you the basics of git.
- [Interactively learn GIT](https://learngitbranching.js.org/)
A git cheat sheet
- [Git cheat sheet](https://education.github.com/git-cheat-sheet-education.pdf)
You can also use git in VSCode (a code editor) to make things easier.
- [VSCode](https://code.visualstudio.com/docs/sourcecontrol/overview)

## Cloning this repository
One of the first things you will want to do is clone this repository.
This essentially copies the repository (code, examples, etc.) from GitHub to you and will allow you to make changes to the code.

To clone this repository, you will need to have an SSH key set up on your GitHub account. Follow these steps to set up an SSH key:
- [Generate an SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) on your machine.
- [Add the SSH key to your GitHub account](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)

Now you can clone the repository, to do this run the following command in your terminal:

```bash
git clone git@github.com:CORDEX-be2/ValEnsPy.git
```

This will create a folder called ValEnsPy in your current directory. You can now navigate to this folder and start working on the code.

## Branching
Branching is a way to work on a feature or bugfix without affecting the main codebase.
Note that ValEnsPy is protected and you will not be able to push directly to the main or dev branch.
Therefore, you are safe to try out branching without worrying about breaking the code.

## Contributing
Now that you have cloned the repository and know some of the basics of git, you can start contributing to the codebase.
To finish setting up your developer environment, follow the instructions in the [installation guide](INSTALL.md).
