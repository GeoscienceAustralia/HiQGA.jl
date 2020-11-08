module SubModule2
    import AbstractModule.foo
    foo(x::Float64) = println("Hello there Float")
end
