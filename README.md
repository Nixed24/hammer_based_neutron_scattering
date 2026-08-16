REQUIRED LIBRARIES:
- numpy
- vmf_deserialiser (included)
- scipy (for constants)


# Hammer-based neutron scattering

This is a somewhat-crazy idea I thought of earlier this year, based off of some uni work I did. The idea is that you use Hammer (```func_brush``` for targets/walls/etc. and ```info_target``` for sources) to build the environment you want to simulate, then the program will take the Hammer VMF and simulate it for you.

For the moment prisms and cubes are the only building blocks you can create your environment with, I plan to make it work with spheres at some point, as most structures you see in the real world can be built by intricate arrangements of those shapes. 

The code contains a technically working example, however there is not yet any absorption mechanism for neutrons as I deemed it more important that all the other components worked well first.

TODO:
- add spheres
- add absorption mechanism for neutrons
- add gamma rays
- add a spatial clamp to iteration alongside the current fixed number of iterations
