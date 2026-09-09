#!/bin/bash

jobs=${NPROC_MAX:-8}
args=( "--clean" "--jobs=${jobs}" "$@" )
if [ "${ASTER_BUILD}" = "debug" ]; then
    args+=( "--timefactor=16.0" )
else
    args+=( "--timefactor=4.0" )
fi


# to be directly readable in the browser
export MESS_EXT="mess.txt"

if [ ! -z "${GITLAB_CI}" ]; then
    echo "+ fetching '${REFREV}' branch..."
    git branch -D ${REFREV} || true
    git fetch --depth=50 origin ${REFREV}
    git branch ${REFREV} FETCH_HEAD
fi
base=$(git merge-base ${REFREV} HEAD)

# check if it only changed testcases
changes=$(git diff --name-status ${base} | grep -v astest/)
if [ -z "${changes}" ]; then
    files=$(git diff --name-status ${base} | awk '{print $2}' | \
        grep astest/ | sed -e 's%astest/%%')
    if [ -z "${files}" ]; then
        printf "no changes detected?!\n"
    else
        ftmp=$(mktemp tmp.list.XXXXXXXX)
        for name in "${files}"; do
            grep -l "${name}" astest/*.export >> ${ftmp}
        done
        flist=$(mktemp tmp.list.XXXXXXXX)
        sort -u ${ftmp} | sed -e 's%astest/%%' -e 's%\.export%%' > ${flist}
        printf "\nchanged testcases:\n"
        cat ${flist}

        # only run these testcases
        args+=( "--testlist=${flist}" )
    fi
fi

if [ "${BUILDTYPE}" = "ci" ] && [ "${CI_JOB_NAME}" != "known_failures_test" ]; then
    args+=( "--exclude-testlist" ".gitlabci/known_failures-pleiade.list" )
fi

run_ctest="./install/bin/run_ctest"
if [ "${OSNAME}" = "win" ]; then
    export LANG=en_EN.UTF-8
    run_ctest="wine ./install/bin/run_ctest.bat"
fi

printf "\nrun_ctest command:\n    "
echo "${run_ctest}"
printf "\nrun_ctest arguments:\n    "
echo "${args[@]}"

printf "\nrunning testcases #1... - $(date)\n"
${run_ctest} "${args[@]}"
iret=$?

if [ ${iret} -ne 0 ]; then
    printf "\nrunning testcases #2 (rerun-failed)... - $(date)\n"
    ${run_ctest} "${args[@]}" --rerun-failed
    iret=$?
fi

if [ ${iret} -ne 0 ]; then
    printf "\nrunning testcases #3 (rerun-failed)... - $(date)\n"
    ${run_ctest} "${args[@]}" --rerun-failed
    iret=$?
fi

if [ -z "${changes}" ] && [ "${BUILDTYPE}" != "nightly" ]; then
    # add mark to skip keywords coverage
    mkdir -p results
    echo "only tests changed" > results/only_tests
fi

# archive gcov data
if [ "${BUILDTYPE}" = "nightly-coverage" ]; then
    printf "\nrunning gcovr... - $(date)\n"
    time gcovr

    printf "\ncreating archive... - $(date)\n"
    tar czf coverage.tgz -C build/mpidebug/debug \
        $(find build/mpidebug/debug -name 'coverage*' -exec basename {} \;)
    mv coverage.tgz results/
fi

printf "\ncreating archives {mess,code}_files.tar.gz... - $(date)\n"
cd results
tar czf mess_files.tar.gz *.${MESS_EXT}
tar czf code_files.tar.gz *.code

printf "\nkeeping output of failed tests in results/failures... - $(date)\n"
mkdir failures
failures=( $(sed -e 's/.*_//g' Testing/Temporary/LastTestsFailed.log 2> /dev/null) )
for test in "${failures[@]}"; do
    mv ${test}.${MESS_EXT} ${test}.code failures/
done
rm -f *.${MESS_EXT} *.code
cd ..

# 'nightly*' runs: archive results files
if [ "${BUILDTYPE:0:7}" = "nightly" ]; then
    cd results

    mc --insecure alias set minio/ ${MINIO_URL} ${MINIO_LOGIN} ${MINIO_PASSWD}
    tdir="${REFREV}"
    [ "${BUILDTYPE}" = "nightly-coverage" ] && tdir="coverage"
    [ "${BUILDTYPE}" = "nightly-sanitize" ] && tdir="sanitize"
    dest=minio/codeaster/devops/ci-${OSNAME}/results/${tdir}/verification
    mc --insecure cp run_testcases.xml ${dest}/
    mc --insecure cp mess_files.tar.gz ${dest}/
    mc --insecure cp code_files.tar.gz ${dest}/
    if [ -f coverage.tgz ]; then
        mc --insecure cp coverage.tgz ${dest}/
        tar xvf coverage.tgz coverage.txt
        today=$(date +%Y-%m-%d)
        mc --insecure cp coverage.txt minio/codeaster/devops/gcov/${today}_coverage.txt
    fi

    cd ..
fi

exit ${iret}
