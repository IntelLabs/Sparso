#!/bin/bash

<<"LICENSE"
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE

# maximal file size to download
max_size=200M

if [[ $# -ge 1 ]]; then
    max_size=$1
fi

max_size_byte=$(numfmt --from=iec ${max_size})

echo "==== Downlaod matrices < `numfmt --to=iec ${max_size_byte}`."

mkdir -p inputs

for line in `cat inputs.list | grep http`; do
    echo
    IFS=',' read -a myarray <<< ${line}
    curr_size_=${myarray[1]}
    curr_size=$(numfmt --from=iec ${curr_size_})
    url=${myarray[0]}
#    curr_size=$(wget "${url}" --spider --server-response -O - 2>&1 | sed -ne '/Content-Length/{s/.*: //;p}')
    target_file="inputs/$(basename $url)"

    download=true
    if ! [[ "$curr_size" =~ ^[0-9]+$ ]]; then
        echo "!!!! error: ${url}" ; continue
    elif [[ -f $target_file ]] && [[ $(stat -c%s "$target_file") == $curr_size ]]; then
        echo "!!!! file exists: ${target_file}"
        download=false
    elif [[ $curr_size -gt $max_size_byte ]]; then
        echo "!!!! skip large file: ${url}, `numfmt --to=iec ${curr_size}`"; continue
    fi

    if [[ $download = true ]]; then
        echo "---- Downloading ${url}, `numfmt --to=iec ${curr_size}`"
        echo
        wget -O $target_file "$url"
    fi

    tar -xf $target_file -C ./inputs/ --strip-components 1
done
