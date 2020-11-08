module SuperModule
import AbstractModule.foo
include("SubModule1.jl")
include("SubModule2.jl")
export foo
end
