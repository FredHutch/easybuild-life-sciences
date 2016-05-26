# HOWTO Build a github pages repo with travis-ci and python-darkslide
---
# Steps (based on 

1. create a repo with a gh-pages branch per [this page] (https://pages.github.com/)
2. create a virtualenv for your project (outside the repo)
3. install python-darkslide: `(your_project)you@yourhost:~/$ pip3 install darkslide`
4. generate a requirements.txt using pip: `(your_project)you@yourhost:~/$ pip3 freeze > requirements.txt`
5. add requirements.txt to your repo
6. create a `.travis.yml` file and add it to your repo - ex:
   !yaml
    language: python
    python:
      - "3.5"
    # command to install dependencies
    install: "pip install -r requirements.txt"
    # command to run tests
    script: darkslide README.md
7. (this is the gnarly part) - on your workstation, run `$ gem install travis` - this may require additional packages like ruby-dev/devel
8. Go get a github personal access token
9. have travis encrypt the token and put in in your .travis.yml file with `travis encrypt GH_TOKEN=<your token> --add`
10. create a regenerate.sh script to update your repo:
   !bash
    #!/bin/bash

    set -o errexit -o nounset

    GIT_USER="bmcgough"
    GIT_EMAIL="bmcgough@fredhutch.org"
    GIT_REPO="FredHutch/easubuild-life-sciences.git"

    if [ "$TRAVIS_BRANCH" != "master" ]
    then
      echo "This commit was made against the $TRAVIS_BRANCH and not master! Not generating."
      exit 0
    fi

    rev=$(git rev-parse --short HEAD)

    # set up a new repo
    git init
    git config user.name "$GIT_USER"
    git config user.email "$GIT_EMAIL"

    # add an upstream remote
    git remote add upstream "https://$GH_TOKEN@github.com/$GIT_REPO"
    git fetch upstream
    git reset upstream/gh-pages

    # touch everything so git will push it
    touch .

    # add, commit, and push
    git add -A .
    git commit -m "regeneration for gh-pages at ${rev}"
    git push -q upstream HEAD:gh-pages

