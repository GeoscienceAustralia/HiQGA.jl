# these are comments
println("Hello!") # you can evaluate with Shift+Enter on this line
#= this is 
also a comment 
=#
println("I'm talking to myself") # try Shift+Enter
## Double hash starts a new code cell in VScode
# you can evaluate it with Alt+Shift+Enter
println("I've evaluated in this cell")
a = ones(3)
println(a)
@info "you can also do" a
println()
println("and also this to show the contents of the last element of a $(a[end]), within a sentence!")
## BIG GOTCHA
b = a # assign alias b to a
b[2] = 2
println("And now b is $b")
println("And now a is $a")
@info b == a
## what you want to do is
c = copy(a)
c[2] = 100
println("a = $a")
println("c = $c")
## plotting
using PyPlot # this is how you call a package
f1 = figure()
plot(rand(10))
f2, ax = plt.subplots(1, 2, sharex=true, figsize=(8,4))
ax[1].imshow(rand(10,10))
ax[2].plot(1:10, 1:10)
# you will use this command a lot
close("all")