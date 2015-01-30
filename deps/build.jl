# See https://github.com/JuliaLang/BinDeps.jl/tree/master/test/testscripts/simplebuild 
# for an example

using BinDeps

const solver_uri = URI("http://to be determined/libpcl_trsolver.tar")
const shlib_ext = BinDeps.shlib_ext

@BinDeps.setup

deps = [
    libpcl_trsolver = library_dependency("libpcl_trsolver", aliases = ["libpcl_trsolver","libpcl_trsolver1","libpcl_trsolver.1"])
]

provides(Sources,solver_uri, libpcl_trsolver,  SHA="to be determined")
provides(BuildProcess,Autotools(libtarget = "libpcl_trsolver.$shlib_ext"),libpcl_trsolver)
@BinDeps.install [:libpcl_trsolver => :jl_libpcl_trsolver]


