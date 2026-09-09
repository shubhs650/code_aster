#!/bin/bash
set -e

TAG=$(git describe --tags --abbrev=0 2>/dev/null || echo "no-tag")
# get the last tag (excluding names)
if [ "${TAG}" != "no-tag" ]; then
    TAG=$(git tag --points-at $(git rev-list -n 1 "${TAG}") | egrep -v '^[a-z]+$' | sort -n | tail -1)
fi
REVISION=$(git rev-parse HEAD)
DISTANCE=$(git rev-list "${TAG}"..HEAD --count 2>/dev/null || echo "0")
DATE=$(git show -s --format=%cd --date=format:%d/%m/%Y HEAD)
BRANCH=$(git branch --show-current || echo "detached")

if [ -z "${BRANCH}" ]; then
    BRANCH="${CI_COMMIT_REF_NAME}"
fi
FROM_BRANCH="${CI_DEFAULT_BRANCH}"

echo "Version info:"
echo "($TAG, $REVISION, $BRANCH, $DATE, $DISTANCE)"

# GitLab dotenv
cat <<EOF > version.env
VERSION_TAG=$TAG
VERSION_REVISION=$REVISION
VERSION_BRANCH=$BRANCH
VERSION_DATE=$DATE
VERSION_DISTANCE=$DISTANCE
EOF

echo "+ setting pkginfo..."
TAGSPLIT=$(echo $TAG | awk -F. '{printf("(%d, %d, %d)\n", $1, $2, $3);}')
cat << EOF > code_aster/pkginfo.py
pkginfo = (${TAGSPLIT}, '${REVISION}', '${BRANCH}', '${DATE}', '${FROM_BRANCH}', ${DISTANCE}, [])
EOF
cat code_aster/pkginfo.py
