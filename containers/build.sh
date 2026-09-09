#!/bin/bash -ex
# because it may use dash...

mkdir -p ${workdir}
cd ${workdir}

echo "+ checking context..."
cd ${workdir}/src
ls -la

echo "+ creating pkginfo..."
if [ ! -f code_aster/pkginfo.py ]; then
    echo "${PKGINFO}" > code_aster/pkginfo.py
else
    echo "already exists (usually from artifacts)"
fi
cat code_aster/pkginfo.py

echo "+ running configure..."
source /opt/public/${DEVTOOLS_COMPUTER_ID}_mpi.sh
./configure --prefix=/opt/codeaster

echo "+ compiling..."
jobs=$(( $(nproc) - 2 ))
make install -j ${jobs}

echo "+ cleaning..."
cd /
rm -rf ${workdir}
