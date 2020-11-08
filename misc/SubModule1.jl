module SubModule1
    import AbstractModule.foo
    foo(x::Int) = println("Hello there Int")
end
