#!/bin/bash

RED='\033[0;31m'
YELLOW='\033[0;33;01m'
RESET='\033[0m' # No Color

set -e

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/hpy"
    exit 1
fi

DIR=$(dirname $0)
HPY=$1

# cd to pypy/module/_hpy_universal/ so we can use relative paths
cd $DIR

# ~~~ helper functions ~~~

indent() {
   sed 's/^/  /'
}

check_dirty() {
    if [[ $(git -C "$HPY" diff --stat) != '' ]]; then
        echo "WARNING! The source hpy repo is dirty"
        echo
    fi
}

check_version_status() {
    # check that the version in hpy/devel/version.py corresponds to the one
    # given by git
    revgit=$(git -C "$HPY" describe --abbrev=7 --dirty --always --tags --long)

    pushd "$HPY/hpy/devel" > /dev/null
    revpy=$(python -c 'import version;print(version.__git_revision__)')
    popd > /dev/null

    if [ "$revgit" != "$revpy" ]
    then
        echo "ERROR: hpy/devel/version.py seems to be outdated"
        echo "  revision reported by git describe: $revgit"
        echo "  revision in hpy/devel/version.py:  $revpy"
        echo
        echo "Please run setup.py build in the hpy repo"
        exit 1
    fi
}

myrsync() {
    rsync --exclude '*~' --exclude '*.pyc' --exclude __pycache__ "$@" 
}

apply_patches() {
    # see also patches/README for more info

    cat > ${DIR}/test/_vendored/conftest.py <<EOF
# AUTOMATICALLY DELETED BY ./update_vendored.sh
EOF
    cat > ${DIR}/test/_vendored/debug/__init__.py <<EOF
# AUTOMATICALLY CREATED BY ./update_vendored.sh
EOF

    fixmes=`ls patches/*FIXME*.patch | wc -l`
    if [ $fixmes -gt 0 ]
    then
        echo -e "${RED}REMINDER: there are ${fixmes} patches marked as FIXME${RESET}:"
        ls -1 patches/*FIXME*.patch | indent
    fi

    for FILE in patches/*.patch
    do
        patch -p4 < $FILE
        if [ $? -ne 0 ]
        then
            echo "${FILE}: patch failed, stopping here"
            echo "See patches/README for more details"
            exit 1
        fi
    done
    echo
}

# ~~~ main code ~~~

check_dirty
check_version_status

myrsync -a --delete ${HPY}/hpy/devel/ ${DIR}/_vendored/hpy/devel/
myrsync -a --delete ${HPY}/test/* ${DIR}/test/_vendored/
apply_patches

echo -e "${YELLOW}GIT status${RESET} of $HPY"
git -C "$HPY" --no-pager log --oneline -n 1
git -C "$HPY" --no-pager diff --stat
echo
echo -e "${YELLOW}HG status${RESET} of pypy"
hg st $DIR
echo
echo -en "${YELLOW}HPy version${RESET}"
cat _vendored/hpy/devel/version.py
