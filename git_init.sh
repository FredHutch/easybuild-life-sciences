#!/bin/bash

set -o errexit -o nounset

BUILD_BRANCH="travis_integration"

if [ "$TRAVIS_BRANCH" != "$BUILD_BRANCH" ]
then
  echo "This commit was made against the $TRAVIS_BRANCH and not $BUILD_BRANCH. Not deploying."
  exit 0
fi

rev=$(git rev-parge --short HEAD)

git init
git config user.name "Ben McGough"
git config user.email "bmcgough@fredhutch.org"

git remote add upstream "https://$GH_TOKEN@github.com/FredHutch/easybuild-life-sciences.git"
get fetch upstream
