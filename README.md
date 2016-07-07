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

# Sparso
*Sparso* is a Julia package for speeding up iterative sparse matrix Julia programs.  
Sparso is part of the High Performance Scripting (HPS) project at Intel Labs.

1. Environment

Sparso is independent of OS, and should work in both Linux and Windows. So far, we
have tested it only on Ubuntu (version >= 14.01, 64 bit desktop).

Julia 0.4.1 has been used for testing.

2. Installation

(1) Getting Julia.

Download a binary from http://julialang.org/downloads/. Untar it and add to PATH.
For example:
	tar xvzf julia-0.4.1-linux-x86_64.tar.gz
	export PATH=path_to_the_julia_directory/bin:$PATH	
	
(2) cd path_to_Julia_package

The path can be found in this way:
	julia   // enter Julia interactive environment
        Pkg.dir()
It shows the path, for example, "/home/username/.julia/v0.4".

(3) git clone https://github.com/IntelLabs/Sparso.git

(4) cd Sparso/scripts; ./setup.sh; cd -

3. Testing correctness

(1) cd Sparso/tests/correctness
(2) julia regression.jl all

You can look at some simple tests there to see how easy it is to use Sparso. Basically,
to accelerate your own iterative sparse matrix applications with Sparso,
add "using Sparso" to your Julia program and add "@acc" in front of calls to functions
that contain the iterative sparse matrix code you'd like to accelerate.

4. Testing performance

Sparso detects sparse matrix operations and redirects where possible to an optimized 
implementation provided by the SpMP library.  Sparso does static and runtime matrix
property discovery in order to select the most performant version of a function on
a matrix having certain properties.  Moreover, Sparso detects constancy of value
and/or structure of matrices allowing different library calls to share an inspector
in the inspector/executor paradigm.  Finally, Sparso does analyses in order to 
determine when to reorder matrices to try to provide higher locality.

(1) cd Sparso/test/perf
(2) julia perf.jl

## Resources

- **GitHub Issues:** <https://github.com/IntelLabs/Sparso.jl/issues>
