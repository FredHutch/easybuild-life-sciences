# HOWTO 
Github with Travis-CI and darkslide
---
# Source
Have Travis-CI auto-generate a slideshow using darkslid and commit it into Github Pages
Steps based on [www.steveklabnik.com](http://www.steveklabnik.com/automatically_update_github_pages_with_travis_example/)
---
# Get darkslide working
   * create a virtualenv for your project outside the repo
   * install python-darkslide: ```(your_project)you@yourhost:~/$ pip3 install darkslide```
   * get your darkslide config working to create the output you want
---
# Set up for Travis-CI
   * generate a requirements.txt using pip: ```(your_project)you@yourhost:~/$ pip3 freeze > requirements.txt```
   * add requirements.txt to your repo
   * create a `.travis.yml` file and add it to your repo
---
# dot(.)travis.yml file
   !YAML
      language: python
      python:
        - "3.5"
      branches:
        only:
        - master
      install:
        - "pip install -r requirements.txt"
        - "bash git_init.sh"
      script:
        - "darkslide README.md -d index.html --relative --copy-theme --theme=default"
        - "darkslide howto.md -d howto.html --relative --copy-theme --theme=default"
      after_success: "bash generate.sh"
      env:
        global:
          secure:
            <encrypted string from travis encrypt goes here>
---
# travis encrypts github token
   * on your workstation, run ```gem install travis``` - this may require additional packages like ruby-dev/devel
   * go get a github personal access token
   * have travis encrypt the token ```travis encrypt GH_TOKEN=<your token> --add``` and add to .travis.yml
   * create a repo with a gh-pages branch per [pages.github.com](https://pages.github.com/)
---
# git_init.sh file
   * create a git_init.sh script to update your repo
   !bash
      #!/bin/bash
      set -o errexit -o nounset
      BUILD_BRANCH="master"
      GIT_USER="bmcgough"
      GIT_EMAIL="bmcgough@fredhutch.org"
      GIT_REPO="FredHutch/easubuild-life-sciences.git"
      if [ "$TRAVIS_BRANCH" != "$BUILD_BRANCH" ]
      then
        echo "commit $TRAVIS_BRANCH and not $BUILD_BRANCH - not generating"
        exit 0
      fi
      # set up a new repo
      git init
      git config user.name "$GIT_USER"
      git config user.email "$GIT_EMAIL"
---
# generate.sh file
   * create a generate.sh script to run darkslide on your markdown:

   !bash
      #!/bin/bash
      set -o errexit -o nounset
      rev=$(git --rev-parse --short HEAD)
      # add, commit, and push
      git add *.html
      git commit -m "regeneration for gh-pages at ${rev}"
      git push -q upstream HEAD:gh-pages

