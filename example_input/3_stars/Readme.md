StarSmasher is also able to handle a 3 star problem. To do so, you have to create first a collision scenario and then se the orbit of the 2 stars that are going to collide. It is not important
how the orbit is.

Then run your simulation and get the first snapshot of it: out0000.sph    , you have to rename it sph.start2u

Now take a relaxed star and name it as sph.start1u


Put them put sph.start1 and 2u inside a folder (like in this example) and compile StarSmasher as we explained in earlier tutorials. What we did? We have first got a collision example and setted it as the system that is going to be hit, and then we toke a single star that is the impactor. In this way, StarSmasher is able to get a 3 star problem.

In this folder we have already procured the stars, in total we have 20000 particles. Compile StarSmasher and try to run it! You can create many scenarios setting sph.input, as you prefer!
