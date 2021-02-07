# math-collocate-perl

Perl module `Math::Collocate` - Collocation with prediction and filtering for scattered data

## Synopsis

```perl
# Include module
use Math::Collocate;

# Create new object using available scattered data
my $this = Math::Collocate->new(sample => \@data);

# Adjust signal and noise variances
$this->sigs($sigs);
$this->sign($sign);

# Apply collocation to create regular grid data
my @result = $this->interpolate({ range => \@range, grid => \@grid });
```

## Requires

+ [List::Util](https://metacpan.org/pod/List::Util)
+ [Math::MatrixReal](https://metacpan.org/pod/Math::MatrixReal)
+ [Moose](https://metacpan.org/pod/Moose)
+ [Scalar::Util](https://metacpan.org/pod/Scalar::Util)
+ [Spreadsheet::Read::Simple](https://github.com/mboljen/spreadsheet-read-simple-perl)

## Installation

To install this module, run the following commands:

```bash
$ perl Makefile.PL
$ make
$ make test
$ make install
```

In order to consider an alternate installation location for scripts,
manpages and libraries, use the `PREFIX` value as a parameter on the
command line:

```bash
$ perl Makefile.PL PREFIX=~/testing
```

## Description

This module implements a least-squares collocation method based on
Moritz [1972].  Collocation is a data processing method which
simultaneously perfoms regression, filtering and spatial prediction based
on given data.  Collocation methods generally assume that measurements
consist of a systematic component (trend surface) and a random component
(signal and noise), i.e.

**w** = **A** **x** + **s** + **n**

where:

- **A** Design matrix from constituent parts of mathematical model
- **x** Vector of unknown parameters (if any)
- **n** Vector of noise components in measurements
- **s** Vector of signal components in measurements
- **w** Vector of measurements

The purpose of least-squares collocation is the prediction and the
separation of the random components **s** and **n** and furthermore,
the estimation of **x** if required.

When the design matrix is known and removed from the measurements,
the task reduces to separate the signal and the noise components:

**z** = **w** - **A** **x** = **s** + **n**.

The constitutents of the vector **z** can be calculated from the
*statistical quadratic form* which is the expression which
weighted least-squares adjustment would minimise:

sqf(**v**) = **v**<sup>T</sup> Cov(**v**,**v**)<sup>-1</sup> **v**.

The vectors **v** contain the distances between the known samples
and the distances to the required point.  The covariance matrix
Cov(**v**,**v**) is the sum of the signal and noise covariance matrices.
The element _ij_ of the *signal* covariance matrix can be computed using
the Gaussian function:

&sigma;<sub>s</sub> ( _s<sub>ij</sub>_ ) = &sigma;<sub>s</sub><sup>2</sup> \* exp( - _s<sub>ij</sub>_<sup>2</sup> / _s_<sub>h</sub><sup>2</sup> )

where _s<sub>ij</sub>_ is the separator factor representing the distance between
the points of interest, and _s_<sub>h</sub> represents some constant control
parameter.  The signal variance &sigma;<sub>s</sub> is set by the user.

If the observed measurements are error-free, then the *noise*
covariance matrix is zero.  Otherwise, a small error is assigned to
each available point _i_:

&sigma;<sub>n</sub> ( _s<sub>ii</sub>_ ) = &sigma;<sub>n</sub><sup>2</sup>

The noise variance &sigma;<sub>n</sub> is set by the user.  By adjusting the
variances &sigma;<sub>s</sub> and &sigma;<sub>n</sub>, the user can control the
separation of the corresponding constituents.

### Export

Nothing.

## Methods

### Constructor

- **new**( _key_ => _value_ )

    The constructor **new** expects a list of key-value pairs and returns the
    reference to the created instance.  The following keys are accepted.

    - **dim** => _int_

        The key **dim** is required if you want to specify the dimension of the problem
        only.  An instance is created without any samples.  Samples can be added later
        using the method **add** described below.

    - **sample** => _arrayref_

        The key **sample** can be used to submit existing data points of arbitrary
        dimension.  If the submitted array reference contains numeric values only,
        the problem is assumed to be 1D and the key **value** is mandatory to submit
        the corresponding values using a reference to an array of the same size.
        If the submitted array reference contains arrays, then the problem is assumed
        to be _N_-dimensional according to the size of the largest array.
        If the key **value** is omitted, the values are taken from the last entry of
        the available arrays and the dimension of the problem is reduced to _N - 1_.

    - **value** => _arrayref_

        The key **value** receives a reference to an array containing numeric values
        that describe the properties to the corresponding samples given by the key
        **sample**.  The user has to take care that the sizes of the arrays referenced
        by **value** and **sample** are the same.

### Object Methods

Here is a list of object methods available.  Object methods are applied to
the object in question, in contrast with class methods which are applied to
a class.

- **size**

    The method **size** returns the number of available samples.

- **sigs**( _value_ )

    The method **sigs** reads or sets the *signal* variance &sigma;<sub>s</sub>
    for the corresponding object.  The default value is `5`.

- **sign**( _value_ )

    The method **sign** reads or sets the *noise* variance &sigma;<sub>n</sub>
    for the corresponding object.  Set the noise variance to zero, if the
    values at the given samples are free of any measurement errors.
    The default value is `1`.

- **add**({ sample => _arrayref_, value => _arrayref_ })

    The method **add** adds a single sample to the corresponding object.  The
    method accepts a single argument, which can either be a hash reference or an
    array reference.

- **prediction**

    The method **prediction** returns the normalized prediction fraction in the
    range `[0..1]`.

- **filter**

    The method **filter** returns the normalized filtering fraction in the
    range `[0..1]`.

- **interpolate**({ range => _arrayref_, grid => _arrayref_ })

    The method **interpolate** applies the collocation algorithm and serves as
    primary method of this class.  The method receives an optional argument which
    can be of any type.  If the method is called in scalar context, the method
    will return a reference to the array containing the results of the
    collocation.  If the method is called in array context, the method will
    return the resulting array directly.

    - **hashref** (see example above)

        If the submitted argument is a hash reference, then the keys **range**
        and **grid** receive the boundaries and the coarseness of the domain of
        interest.

        - **range** => _arrayref_

            The optional key **range** receives an array reference either to scalars
            indicating the lower and upper boundary of the requested range separated
            by a colon (`:`) or to arrays with 2 values, i.e. the lower and the
            upper boudnary itselves.  If several array elements are undefined or the
            array is not defined at all, the range is determined automatically based
            on the minimum and the maximum value given by the samples.

        - **grid** => _int_ | _arrayref_

            The optional key **grid** may receive a positive integer which is applied
            to all dimensions available or it may receive a array of integers which
            indicate the requested grid size for each dimension.  If key is omitted,
            the default grid of `20` is applied.

    - **arrayref**

        If the submitted value is an array reference, the value is assumed to be a
        single sample with size corresponding to the dimension of the object.
        The method will return the resulting value in scalar context.

    - **scalar**

        A scalar value is allowed for objects with dimension 1 only, otherwise an
        error is raised.  The method will return the resulting value as a scalar.

- **newline**( _index_ )

    The method **newline** returns `true` or `false` whether the sample
    referred to by the given _index_ closes a block, i.e. is a multiple of
    the grid points along the first dimension.  This method is particularly
    useful to insert a newline when the output data shall be forwarded
    to **gnuplot**.

## Subroutines

The following routines are provided.

- **isequal**( _value1_, _value2_ [, _tol_ ] )

    This subroutine compares two floating-point numbers _value1_ and
    _value2_ and checks whether they are equal within a reasonable
    tolerance _tol_.

## References

+ Höpcke, W. (1980):

  *Fehlerlehre und Ausgleichsrechnung*.  De Gruyter, Pages 211 f.

- Moritz, H. (1972):

  *Advanced least-squares methods*.  Report No. 175, Department of Geodetic Science, Ohio State University.  132 Pages.

- Ruffhead, A. (1987):

  *An introduction to least-squares collocation*.  Survey Review, Volume 29, Number 224, Pages 85-94.

## Acknowledgements

My father Joachim Boljen for guiding me the way through the mists of statistical geodesy for filtering, prediction and modeling.

## Copyright and Licnse

MIT License

Copyright (c) 2020, 2021 Matthias Boljen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
