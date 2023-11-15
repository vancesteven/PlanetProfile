# Contributing to PlanetProfile
PlanetProfile is open source software. See the [LICENSE](https://github.com/vancesteven/PlanetProfile/blob/main/LICENSE) file and read the [README.md](https://github.com/vancesteven/PlanetProfile/blob/main/README.md) file first. Please also review our community [CODE_OF_CONDUCT.md](https://github.com/vancesteven/PlanetProfile/blob/main/CODE_OF_CONDUCT.md).

We welcome input from everyone! If you are interested to contribute to PlanetProfile, please read this document carefully. If you have any questions, please contact the software Maintainer at the email listed below. PlanetProfile leadership all maintain this software on a part-time, volunteer basis. We endeavor to respond to any inquiry within 3 business days of receipt. If you don't receive a response by then, feel free to ping again.

## Governance
* Owner: Dr. Steven D. Vance - steven.d.vance@jpl.nasa.gov
* Maintainer: Dr. Marshall J. Styczinski - itsmoosh@gmail.com
* Lead developers: Dr. Steven D. Vance - steven.d.vance@jpl.nasa.gov, Dr. Marshall J. Styczinski - itsmoosh@gmail.com, Dr. Mohit Melwani Daswani - mohit.melwani.daswani@jpl.nasa.gov

## Contribution guidelines
PlanetProfile is designed for use in scientific applications. As such, we hold the software to a high standard for usability, compatibility, and accuracy. We review all contributions to ensure the quality of all improvements to the software. Proposed changes may be rejected if they cause runtime errors, physical inaccuracy, etc. In the event a change is rejected, we will work with the author to find a solution. Command line operations are assumed in the instruction below, but you can substitute another method if it has the same effect.

We have a Slack workspace that we use for development discussions. You will be invited to join the Slack once you express interest in contributing.

### Opening issues
This software is in active development. For any issues you find, please use [the GitHub issue reporting tool](https://github.com/vancesteven/PlanetProfile/issues) to document the problem and any solution you may be aware of. Including filenames and line numbers makes issues much easier to track, so please identify these whenever possible. To make sure your issue is addressed promptly, follow the checklist in the issue template. The template appears when you attempt to submit an issue.

### Making commits
As PlanetProfile is written at least partially in Matlab, it's up to us contributors to test that our commits function properly.

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

## For developers: Packaging
The python version of PlanetProfile is listed on the Python Package Index, PyPI. Updating the package and release version require access permissions (which can be granted by Dr. Vance) and the following steps:
1. Update documentation
   1. Navigate to the docs/ folder. From the top-level directory: `cd docs`
   1. Generate HTML files: `rm -rf stubs/ && make clean && make html`. Make a note of any warnings.
   1. Open `docs/_build/html/index.html` in a browser and confirm that everything looks as it should.
   1. Resolve any warnings or issues with the documentation web pages before committing changes.
1. Create pypi package
   1. Delete previous package files from dist/ folder.
   1. Update version number in both `setup.py` and `PlanetProfile/Utilities/PPverNum.txt`
   1. Construct the package from `setup.py` and `MANIFEST.in` instructions: `python -m build`
1. Upload pypi package and GitHub/Zenodo release
   1. Upload to PyPI: `python -m twine upload dist/* --verbose`
      1. For Username, enter verbatim: `__token__` and hit enter.
      1. For Pass, copy to clipboard and paste into the terminal your PyPI API token and hit enter. The token is what you saved during your PyPI account setup.
   1. Upload GitHub/Zenodo release:
      1. Navigate to <https://github.com/vancesteven/PlanetProfile/releases/new>
      1. Drag-and-drop or upload the new `dist/PlanetProfile-#.#.#.tar.gz` file
      1. Enter the new version number, preceded by "v", in the "Choose a tag" menu at the top. Target should be "main".
      1. Give the release a title that briefly indicates what's new in this version.
      1. Add a couple sentences describing the changes in slightly more detail in the "Describe this release" box. See previous releases for examples.
1. Update the mirror in the NASA-planetary-science org
   1. Navigate to the mirror repo local directory (e.g. `cd ~/NASApp/`)
   1. Update the mirror: `git remote update && git push --mirror`

### Initial setup
Going through these steps for the first time requires some initial setup. Follow the steps below the first time you are ready to create a new release. You will need write permissions for:
1. The main repo at <https://github.com/vancesteven/PlanetProfile>
2. The PyPI package at <https://pypi.org/project/PlanetProfile/>
3. The mirror repo at <https://github.com/NASA-Planetary-Science/PlanetProfile>

Contact project personnel to get write permissions if required.

#### Creating and uploading the PyPI package
1. Create an account on PyPI: <https://pypi.org/account/register/>
   1. Confirm your account by email
   1. Generate and save (in a safe place!) your upload authentication token. You will need it every time you update the package on PyPI.
1. Ensure pip is the latest version and that build and twine are installed: `pip install --upgrade pip build twine`

#### Testing the documentation build
Install prerequisites with pip: `pip install sphinx sphinxcontrib-matlabdomain sphinxcontrib.apidoc sphinx-rtd-theme myst-parser`
* Sphinx is a tool for reading docstrings in the code and converting them to documentation pages.
* sphinxcontrib-matlabdomain provides support for Matlab docstrings.
* sphinx-rtd-theme provides the 'sphinx_rtd_theme' option (note the underscores) for the Read-the-Docs theme
* myst-parser provides a conversion tool for parsing README files written in Markdown (.md files)

#### Preparing to update the mirror repository
1. Navigate to where you want to local GitHub files to live (e.g. your home directory)
2. Clone a mirror of the main repo: `git clone --mirror git@github.com:vancesteven/PlanetProfile.git`
3. Rename the mirror to avoid conflicts with the main repo. Recommended: `mv PlanetProfile.git NASApp`
4. Navigate to the new folder, e.g. `cd NASApp`
5. Edit the config file using `git config -e` to point refs from the main repo to the mirror. That is:\
   Change the line\
      `fetch = +refs/*:refs/*`\
   to instead be\
      `fetch = +refs/tags/*:refs/tags/*`\
      `fetch = +refs/heads/*:refs/heads/*`\
      `push = +refs/tags/*:refs/tags/*`\
      `push = +refs/heads/*:refs/heads/*`\
      `pushurl = git@github.com:NASA-Planetary-Science/PlanetProfile.git`\
   Then save and exit.

### Building the docs from scratch
In the event of a major restructure, Sphinx is not well set up to reformat the docs without a lot of manual intervention. Here are the steps to redraw the documentation pages from scratch.
1. Install sphinx prerequisites and extensions [listed above](#testing-the-documentation-build).
1. Copy any custom or edited documentation files, including: 
   1. `conf.py`
   1. `index.rst`
   1. `readme.rst`
   1. `requirements.txt`
   1. `pythonModules.txt`
   1. `.gitignore`
   1. `.nojekyll`
   1. `misc/` and its contents
   1. `SPICE/` and its contents
1. Delete the entire `docs` directory
1. Navigate to the repo top-level directory and run the following command:\
   `sphinx-quickstart docs`
   1. Accept the default and do not separate source and build directories
   1. Enter `PlanetProfile` for the project name
   1. Enter `Vance` for the author name (we will overwrite this with conf.py in a moment, but we have to write something)
   1. Accept the defaults for version number and language
1. Delete the unnecessary `docs/_templates` directory
1. Copy the custom files you saved back into the `docs/` directory, replacing the existing files
1. Navigate into the `docs` directory with `cd docs` and generate pages with `make html`

Subsequent builds need only run `rm -rf stubs/ && make clean && make html` from the `docs/` directory.
