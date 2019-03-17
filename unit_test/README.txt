To run unit tests from MATLAB:
  

Ensure that SMTB is on your path, then

```matlab
>> RunAllTests
```

will also generate a [Cobertura](http://cobertura.github.io/cobertura/) format coverage report which is viewable using services like [codecov](https://codecov.io/gh/petercorke/spatial-math).

This file is invoked from Travis CI via the `matlab-runner`.
