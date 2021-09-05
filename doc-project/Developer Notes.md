# Neuropixel Utils Developer Notes

## Repo layout

* `Mcode/` – Matlab source code for the Neuropixel Utils library itself
* `doc-project/` – developer-oriented doco for people working on the library itself
* `doc-src/` – Source for user-facing doco
* `docs/` – Generated user-facing doco and GitHub Pages site
  * This is generated from `doc-src/`. Don't edit anything here directly!
* `dev-kit/` – Developer tools for working on the library
* `examples/` – Example code to go with the library, and for building into the `docs/` stuff
* `map_files/` – ??? Are these example data sets or part of the library itself?
* Ephemeral directories (not checked in to Git)
  * `build/` – Where intermediate munged source code goes as part of the build
  * `dist/` – Final distribution artifacts appear here

## Working with the repo

When working on Neuropixel Utils as a developer, you'll want both the NeuroPixel Utils library and its dev kit loaded on to your path. Run this to load it up:

```matlab
% Edit this to reflect where you put your copy of the repo
npxutilsRepoPath = '~/local/repos/neuropixel-utils';

addpath(fullfile(npxutilsRepoPath, 'dev-kit'))
load_npxutils
```

It can be handy to stick that code in a Favorite in your Matlab. Select Home > Favorites > New Favorite from the toolstrip.

## Workflow notes

### ChangeLog

The [ChangeLog](https://keepachangelog.com/) for Neuropixels Utils is at `CHANGES.md` in the root of the repo.

Whenever you make a significant, user-visible change to the library (especially if it's an API change!), make a note about it in `CHANGES.md` right away. Don't wait until release time to do this; it's hard to remember 

## Building the project

### Building the web site

The project web site is a GitHub Pages site, but it uses MkDocs instead of the default Jekyll as its static site generator.

To build the site with MkDocs:

First, get MkDocs installed. With Homebrew:

```bash
brew install mkdocs
```

Or generically, with Anaconda Python:

```bash
conda create --name mkdocs pip
conda activate mkdocs
pip install mkdocs mkdocs-material
```

Then, to build the docs using MkDocs:

```bash
mkdocs serve # test local
mkdocs gh-deploy # deploy to master
```

### Building the distribution

With the dev-kit loaded as described above, run this in Matlab to build the library as a Matlab Toolbox `.mltbx` file:

```matlab
npxutils_make toolbox
```

And run this to build the distribution zip files:

```matlab
npxutils_make dist
```

This will produce the distribution Matlab Toolbox .mltbx file.

## Style Guide

### Prose Style

* "Matlab", not "MATLAB"
  * Yes, we know "MATLAB" is the official styling; we just don't like it.
* Oxford commas
* Full stops (periods) at the end of list items that are sentences or similar

### Code Style

* 4-space indents
* Auto-format everything using the Matlab Editor
* Naming
  * `camelCase` variable, function, and method names
  * `TitleCase` class names
  * `lowercase` package names
  * Method dispatch objects are named `this`
* Markdown
  * MkDocs-flavored Markdown
    * (Andrew isn't exactly sure what this is, actually!)
  * Try to be Markdownlint-clean
