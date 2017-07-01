# Cubic Spline

This project was made as an extension of my homework task at Numerical analysis course at schoolyear 2015/16.

Project consists of `CubicSpline` class which interpolates function using [cubic spline interplation](http://mathworld.wolfram.com/CubicSpline.html).

### Compiling

You need to have installed C++11 compiler or later version to compile this project.
When you have it, just run make from root folder for compiling.

### Usage

Just run compiled executable and follow given instructions.

Alternatively, run executable with given path to file with spline data (example is `spline_data.txt`). File format should be following:
* first `n` lines consisting of two numerical values separated by space represent `n` datapoints which will be interpolated
* one empty line
* line with two numerical values - second derivatives of function in it's first and last node, respectively.

### Development

This is more or less done project. If you feel it should be expanded, feel free to push an update.

### Author
Tomislav Bukić


License
----

This project is under MIT license.
