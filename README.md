[![Build Status](https://travis-ci.com/petercorke/spatial-math.svg?branch=master)](https://travis-ci.com/petercorke/spatial-math)
![Coverage](https://codecov.io/gh/petercorke/spatial-math/branch/master/graph/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/petercorke/robotics-toolbox-matlab/graphs/commit-activity)
[![GitHub stars](https://img.shields.io/github/stars/petercorke/spatial-math.svg?style=social&label=Star&maxAge=2592000)](https://GitHub.com/petercorke/robotics-toolbox-matlab/stargazers/)
 
# Spatial Math Toolbox for MATLAB&reg;

## Synopsis

This Toolbox contains functions and classes to represent orientation and pose in 2D and 3D (SO(2), SE(2), SO(3), SE(3)) as orthogonal and homogeneous transformation matrices, quaternions, twists, triple angles, and matrix exponentials. The Toolbox also provides functions for manipulating these datatypes, converting between them, composing them, transforming points and graphically displaying them.

Much of this Toolbox functionality was previously in the [Robotics Toolbox for MATLAB](https://github.com/petercorke/robotics-toolbox-matlab).

Advantages of the Toolbox are that:

  * the code is mature and provides a point of comparison for other implementations of the same algorithms;
  * the routines are generally written in a straightforward manner which allows for easy understanding, perhaps at the expense of computational efficiency. If you feel strongly about computational efficiency then you can always rewrite the function to be more efficient, compile the M-file using the MATLAB compiler, or create a MEX version;
  * source code is available for the benefit for understanding and teaching.
  
Comprehensive detail in the [PDF-format manual (193 pages)](http://www.petercorke.com/SMTB/spatialmath.pdf).

## Code Example

```matlab
>> R = rotx(0.2)  % SO(3) rotation matrix
R =
    1.0000         0         0
         0    0.9801   -0.1987
         0    0.1987    0.9801
```

which we could animate simply as
```matlab
>> tranimate(R)
```

![animation from tranimate()](doc/figs/tranimate.gif)

## What's new

* Continuous intergration using [Travis CI](travis-ci.com) and [codecov](codecov.io)
* Support for spatial vector notation (Featherstone's 6D vectors)
* `prod()` method for all `RTBPose` subclasses and `Twist`
* Factored out of the [Robotics Toolbox for MATLAB](https://github.com/petercorke/robotics-toolbox-matlab).  RTB now contains only robotics specific code, eg. for manipulator arms and mobile robots.

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
