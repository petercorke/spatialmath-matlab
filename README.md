 #Spatial Math Toolbox for MATLAB&reg;.

## Synopsis

This Toolbox contains functions and classes to represent orientation and pose in 2D and 3D (SO(2), SE(2), SO(3), SE(3)) as matrices, quaternions, twists, triple angles, and matrix exponentials. The Toolbox also provides functions for manipulating and converting between datatypes such as vectors, homogeneous transformations and unit-quaternions which are necessary to represent 3-dimensional position and orientation.

This Toolbox has been factored out of the [Robotics Toolbox for MATLAB](https://github.com/petercorke/robotics-toolbox-matlab).

Advantages of the Toolbox are that:
  * the code is mature and provides a point of comparison for other implementations of the same algorithms;
  * the routines are generally written in a straightforward manner which allows for easy understanding, perhaps at the expense of computational efficiency. If you feel strongly about computational efficiency then you can always rewrite the function to be more efficient, compile the M-file using the MATLAB compiler, or create a MEX version;
  * source code is available for the benefit for understanding and teaching.

## Code Example

```matlab
>> rotx(0.2)  % SO(3) rotation matrix
ans =
    1.0000         0         0
         0    0.9801   -0.1987
         0    0.1987    0.9801
```

## What's new

* Support for spatial vector notation (Featherstone's 6D vectors)
* `prod()` method for all `RTBPose` subclasses and `Twist`

## Online resources:

* [Home page](http://www.petercorke.com)
* [Robotics, Vision & Control (book)](http://petercorke.com/wordpress/rvc)
* [Discussion group](http://groups.google.com/group/robotics-tool-box?hl=en)

Please email bug reports, comments or code contribtions to me at rvc@petercorke.com

## Octave

The functions, but not the classes, should work fine with Octave 4.x.

## Contributors

Contributions welcome.  There's a user forum at [http://tiny.cc/rvcforum](http://tiny.cc/rvcforum) for this Toolbox and also
[Robotics Toolbox for MATLAB](https://github.com/petercorke/robotics-toolbox-matlab).

## License

This toolbox is released under MIT Licence.
