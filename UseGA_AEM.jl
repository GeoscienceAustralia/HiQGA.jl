module UseGA_AEM
using Libdl
using Cxx 

mutable struct emforward
    TLM  :: Cxx.CxxCore.CppPtr
    THM  :: Cxx.CxxCore.CppPtr
    PX   :: Float64
    PY   :: Float64
    PZ   :: Float64
    SXLM :: Array{Float64, 1}
    SYLM :: Array{Float64, 1}
    SZLM :: Array{Float64, 1}
    SXHM :: Array{Float64, 1}
    SYHM :: Array{Float64, 1}
    SZHM :: Array{Float64, 1}
end    

function init_GA_AEM()
    #CXXJL_ROOTDIR should be set at least when first time installing Cxx
    #export CXXJL_ROOTDIR=/apps/gcc/5.2.0
    
    # Importing shared library and header file
    path_to_lib = "/g/data1a/zk34/ga-aem/julia"
    addHeaderDir(path_to_lib, kind=C_System)

    ga_aem_path = "/home/547/ar0754/zk34/ga-aem"
    
    # Julia compiler also needs to find headers in these dirs
    addHeaderDir("/g/data1a/zk34/ga-aem/src", kind=C_System)
    addHeaderDir("/g/data1a/zk34/ga-aem/submodules/cpp-utils/src/", kind=C_System)
    addHeaderDir("/apps/fftw3/3.3.7-gcc/include", kind=C_System)
    
    cxxinclude("gatdaem1d_julia.h")
    Libdl.dlopen(path_to_lib * "/gatdaem1d_julia.so", Libdl.RTLD_GLOBAL)
    
    # Create the cTDEmSystem class object
    stmfile = "/g/data1a/zk34/ga-aem/examples/bhmar-skytem/stmfiles/Skytem-LM.stm"
    TLM = @cxxnew cTDEmSystem(pointer(stmfile));
    stmfile = "/g/data1a/zk34/ga-aem/examples/bhmar-skytem/stmfiles/Skytem-HM.stm"
    THM = @cxxnew cTDEmSystem(pointer(stmfile));
    
    # Allocate reusable storage for outputs
    PX = 0.0
    PY = 0.0
    PZ = 0.0
    nwindowsLM = @cxx TLM->NumberOfWindows;
    nwindowsHM = @cxx THM->NumberOfWindows;
    SXLM = Array{Float64, 1}(undef,nwindowsLM)
    SYLM = Array{Float64, 1}(undef,nwindowsLM)
    SZLM = Array{Float64, 1}(undef,nwindowsLM)

    SXHM = Array{Float64, 1}(undef,nwindowsHM)
    SYHM = Array{Float64, 1}(undef,nwindowsHM)
    SZHM = Array{Float64, 1}(undef,nwindowsHM)
    
    # call it once bug
    # Set the input arrays
    geometryLM   = [35.0,0.0,0.0,0.0,-17.0,0.0,+2.0,0.0,0.0,0.0]
    geometryHM   = [35.0,0.0,0.0,0.0,-17.0,0.0,+0.2,0.0,0.0,0.0]
    conductivity = [0.01, 1.0, 0.001]
    thickness    = [20.0, 40.0]
    nlayers      = length(conductivity)

    @cxx TLM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),
                           pointer(geometryLM),Ref(PX),Ref(PY),Ref(PZ),
                           pointer(SXLM),pointer(SYLM),pointer(SZLM))
    @cxx THM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),
                           pointer(geometryHM),Ref(PX),Ref(PY),Ref(PZ),
                           pointer(SXHM),pointer(SYHM),pointer(SZHM))


    emforward(TLM, THM, PX, PY, PZ, SXLM, SYLM, SZLM, SXHM, SYHM, SZHM)
end

function(em::emforward)()
    
    # Set the input arrays
    geometryLM   = [35.0,0.0,0.0,0.0,-17.0,0.0,+2.0,0.0,0.0,0.0]
    geometryHM   = [35.0,0.0,0.0,0.0,-17.0,0.0,+0.2,0.0,0.0,0.0]
    conductivity = [0.01, 1.0, 0.001]
    thickness    = [20.0, 40.0]
    nlayers      = length(conductivity)

    @cxx em.TLM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),
                           pointer(geometryLM),Ref(em.PX),Ref(em.PY),Ref(em.PZ),
                           pointer(em.SXLM),pointer(em.SYLM),pointer(em.SZLM))
    @cxx em.THM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),
                           pointer(geometryHM),Ref(em.PX),Ref(em.PY),Ref(em.PZ),
                           pointer(em.SXHM),pointer(em.SYHM),pointer(em.SZHM))
end

end
