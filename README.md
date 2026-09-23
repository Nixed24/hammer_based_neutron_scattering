REQUIRED LIBRARIES:
- numpy
- vmf_deserialiser (included)
- scipy (for constants)
REQUIRED APPLICATIONS:
- Hammer World Editor (should be bundled with any Source Engine game's Windows binaries)

# Hammer-based neutron scattering

This is a somewhat-crazy idea I thought of earlier this year, based off of some uni work I did. The idea is that you use Hammer (```func_brush``` for targets/walls/etc. and ```info_target``` for sources) to build the environment you want to simulate, then the program will take the Hammer VMF and simulate it for you.

For the moment prisms and cubes are the only building blocks you can create your environment with, I plan to make it work with spheres at some point, as most structures you see in the real world can be built by intricate arrangements of those shapes. 

The code contains a technically working example, however there is not yet any absorption mechanism for neutrons as I deemed it more important that all the other components worked well first.

HOW TO USE:

The 5 currently supported material-texture pairs are:

Material: Water | Texture: liquids/water_pretty1
Material: Air | Texture: None, will be the default material assumed if a particle isn't in a solid (`material_fallback` in the program file)
Material: Lead | Texture: building_template/building_template021a
Material: Graphite | Texture: nature/dirtfloor003a
Material: Anti | Texture: tools/toolsnodraw // Not working yet, will be added next update. For now, use the `is_anti_solid` key value.
"Anti" is a special material that makes solids with this texture act as a subtraction to any other solid that interpenetrates it, see https://en.wikipedia.org/wiki/Constructive_solid_geometry for more details

By default, solid regions negated using the Anti material will be treated as Air, this is controlled by the `material_anti` variable.

In order to create an environment to simulate, open Hammer and construct it using its tools out of the supported textures, or any you have added yourself.
If you haven't used Hammer before, try following the guide at https://developer.valvesoftware.com/wiki/Getting_Started.

Make sure all solids that you make are `func_brush` entities that have only one solid.

In order to make a solid a detector, you must turn off SmartEdit (found near the top and slight right of the entity keyvalue window) in the keyvalue editor. Then, add a keyvalue with key name `is_detector` with the key value set to `yes`.

In order to make a solid a subtractor solid (see Material: Anti), add a keyvalue with key name `is_anti_solid` and set it to `yes`.

In order to make a point source, create an `info_target` point_entity and name it according to the syntax {particle_type}_nickname, e.g `neutron_test_source` for a neutron source.
Then, add keyvalue pairs with key names `amount` and `initial_particle_speed`, to specify the amount of particles for the point source to simultaneously project and the speed to project them at.
Specify `initial_particle_speed` values in units of metres per second.

Use the `angles` keyvalue to specify the direction for the source to project particles in.

Once you have created your environment in Hammer, save it as a VMF, and then change the `VMF_FILENAME` variable to point to your saved VMF.

Then, run the program, and the results of the simulation should be printed in the interpreter.


TODO:
- add spheres
- add absorption mechanism for neutrons
- add gamma rays
- add a spatial clamp to iteration alongside the current fixed number of iterations
