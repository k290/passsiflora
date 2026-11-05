#!/bin/bash
# Fetches the upstream source Passiflora is built from (docs/BUILDING.md, step 1): drivere/get-decics-numberfields at commit
# 638ca4e, only GetDecicsSrc, Test (the work units and their truth files) and PARI 2.8-1711.
#   fetch-upstream.sh [directory]      (default: upstream)
# The partial clone can break off on a slow line ("early EOF"); the checkout is then retried over HTTP/1.1.
set -e
DIR=${1:-upstream}
git clone -q --filter=blob:none --no-checkout https://github.com/drivere/get-decics-numberfields "$DIR"
cd "$DIR"
git sparse-checkout init --cone
git sparse-checkout set GetDecicsSrc Test Pari/pari-2.8-1711-ge5c317c
git checkout -q 638ca4e || git -c http.version=HTTP/1.1 checkout -q -f 638ca4e
echo "upstream ready: $DIR ($(git rev-parse --short HEAD))"
