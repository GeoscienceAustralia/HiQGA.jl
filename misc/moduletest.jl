any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using Revise, SuperModule
