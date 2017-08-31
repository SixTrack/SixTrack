#! /usr/bin/env bash
# Small script to make adding tags a bit less annoying
# This is meant to be ran in a clone linked to the main SixTrack repository,
# by someone who has push access directly there.
#
# This script will create the tag and date it correctly;
# it is then up to the user to create a release from this tag.
#
# K.Sjobak, 31/08/2017
#

if [ "$#" -ne 2 ]; then
    echo "Usage:"
    echo "./$1 v1.2.3 HASH"
    exit
fi

tag=$1
hash=$2

echo "Creating tag='$tag'"
echo "pointing to commit='$hash'"
echo

git tag -a $tag $hash -m "Version $tag"

COMMIT_HASH=$(git rev-list -1 $tag)
GIT_COMMITTER_DATE="$(git show $COMMIT_HASH --format=%aD | head -1)" git tag -a -f $tag -m"$tag" $COMMIT_HASH

git push --follow-tags
