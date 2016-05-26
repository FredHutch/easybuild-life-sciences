#!/bin/bash

set -o errexit -o nounset
set -x

# revision
rev=$(git rev-parse --short HEAD)

# switch to gh-pages
git reset upstream/gh-pages

# touch everything so git will push it
touch .

# add, commit, and push
git add -A .
git commit -m "regeneration for gh-pages at ${rev}"
git push -q upstream HEAD:gh-pages
