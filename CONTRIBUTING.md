# Contributing to PlanetProfile
PlanetProfile is open source software. See the [LICENSE](https://github.com/vancesteven/PlanetProfile/blob/main/LICENSE) file and read the [README.md](https://github.com/vancesteven/PlanetProfile/blob/main/README.md) file first. Please also review our community [CODE_OF_CONDUCT.md](https://github.com/vancesteven/PlanetProfile/blob/main/CODE_OF_CONDUCT.md).

We welcome input from everyone! If you are interested to contribute to PlanetProfile, please read this document carefully. If you have any questions, please contact the software Maintainer at the email listed below. PlanetProfile leadership all maintain this software on a part-time, volunteer basis. We endeavor to respond to any inquiry within 3 business days of receipt. If you don't receive a response by then, feel free to ping again.

## Governance
* Owner: Dr. Steven D. Vance - steven.d.vance@jpl.nasa.gov
* Maintainer: Dr. Marshall J. Styczinski - marshall.j.styczinski@jpl.nasa.gov
* Lead developers: Dr. Steven D. Vance - steven.d.vance@jpl.nasa.gov, Dr. Marshall J. Styczinski - marshall.j.styczinski@jpl.nasa.gov, Dr. Mohit Melwani Daswani - mohit.melwani.daswani@jpl.nasa.gov

## Contribution guidelines
PlanetProfile is designed for use in scientific applications. As such, we hold the software to a high standard for usability, compatibility, and accuracy. We review all contributions to ensure the quality of all improvements to the software. Proposed changes may be rejected if they cause runtime errors, physical inaccuracy, etc. In the event a change is rejected, we will work with the author to find a solution. Command line operations are assumed in the instruction below, but you can substitute another method if it has the same effect.

We have a Slack workspace that we use for development discussions. You will be invited to join the Slack once you express interest in contributing.

### Opening issues
This software is in active development. For any issues you find, please use [the GitHub issue reporting tool](https://github.com/vancesteven/PlanetProfile/issues) to document the problem and any solution you may be aware of. Including filenames and line numbers makes issues much easier to track, so please identify these whenever possible. To make sure your issue is addressed promptly, follow the checklist in the issue template. The template appears when you attempt to submit an issue.

### Making commits
As PlanetProfile is currently written in Matlab, it's up to us contributors to test that our commits function properly.

Before *any and all* commits or pull requests to Matlab code, change all CALC_NEW flags in config.m to 1, then run the test function PPtest in Matlab. This function must run without errors for your changes to be accepted.

If you would like to be considered for direct commit permissions, please contact Dr. Vance at the email listed above.

## Commit checklist
If you would like to propose changes to PlanetProfile, please use one of the following procedures.

### New contributors
Please read [GitHub docs about pull requests](https://docs.github.com/en/free-pro-team@latest/github/collaborating-with-issues-and-pull-requests/about-pull-requests) before continuing.
1. Fork the main repository to your own GitHub profile.
   1. On the [main repository page](https://github.com/vancesteven/PlanetProfile), click the "Fork" button (usually at top-right).
1. In your forked version, create a new branch with a name appropriate to the changes you wish to make.
1. Download the forked version and apply (don't push yet) only the changes you wish to commit to the main repository.
1. Change all CALC_NEW flags in config.m to 1, then run PPtest in Matlab. Fix any errors before continuing.
1. Push your changes *to the new branch* in your forked version.
1. On the web page for your forked version (at `https://github.com/`yourUsername`/PlanetProfile`), click the "Pull Request" link just next to the colored "Code" button.
1. Write a descriptive commit message detailing the changes.
   1. A detailed comment is encouraged but not required. For example, describing the reasoning and content of your changes will make it easier for a maintainer to verify your contributions.
1. A maintainer will review your commit and inform you through GitHub if further action is required. Thank you for contributing!
### Contributors with commit permissions
1. Align with the repository.
   1. Save your local edits with `git stash`
   1. Ensure your local version of the repository is up-to-date with `git pull`
   1. Re-apply your local edits with `git stash pop`
   1. Use `git status` to identify any files with merge conflicts.
   1. Resolve merge conflicts by opening affected files and searching for `====`, which is used to separate the conflicting contributions.
1. Test your edits for runtime errors.
   1. Change all CALC_NEW flags in config.m to 1, then run PPtest in Matlab. Fix any errors before continuing.
1. Apply your changes to the repository.
   1. Use `git add <filepath1> <filepath2> ...` to mark your desired changes for commit.
   1. Again use `git status` to check that everything is as you intended it.
   1. Use `git commit -m "Commit message"` to ready your changes for application.
   1. Use `git push` to apply your changes.
