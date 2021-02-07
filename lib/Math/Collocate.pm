# Copyright 2020 Matthias Boljen.  All rights reserved.
#
# Created:         Mi 2020-10-28 07:51:16 CET
# Last modified:   Di 2021-02-02 17:54:42 CET
#
# This program is free software; you can redistribute it and/or
# modify it under the same terms as Perl itself.

package Math::Collocate;

use strict;
use warnings;
use 5.010;

require Exporter;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@ISA = qw(Exporter);
@EXPORT = qw();
%EXPORT_TAGS = ( 'all' => [ qw() ] );
@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

$VERSION = 0.01;

use Carp;
use List::Util qw(all any max min sum uniqnum);
use Math::MatrixReal;
use Moose;
use Regexp::Common qw(number);
use Scalar::Util qw(reftype);

has 'sample' => (
    isa => 'ArrayRef[Num]|ArrayRef[ArrayRef[Num]]',
    is  => 'ro',
);

has 'value' => (
    isa => 'ArrayRef[Num]',
    is  => 'ro',
);

has 'dim' => (
    isa => 'Int',
    is  => 'ro',
);

has 'grid' => (
    isa => 'ArrayRef[Int]',
    is  => 'ro',
);

has 'sign' => (
    isa     => 'Num',
    is      => 'rw',
    default => 1,
);

has 'sigs' => (
    isa     => 'Num',
    is      => 'rw',
    default => 5,
);

has 'sh' => (
    isa => 'Maybe[Num]',
    is  => 'rw',
);

=encoding utf8

=head1 NAME

Math::Collocate - Collocation with prediction and filtering for scattered data

=head1 SYNOPSIS

    # Include module
    use Math::Collocate;

    # Create new object using available scattered data
    my $this = Math::Collocate->new(sample => \@data);

    # Adjust signal and noise variances
    $this->sigs($sigs);
    $this->sign($sign);

    # Apply collocation to create regular grid data
    my @result = $this->interpolate({ range => \@range, grid => \@grid });

=head1 REQUIRES

List::Util, Math::MatrixReal, Moose, Scalar::Util

=head1 DESCRIPTION

This module implements a least-squares collocation method based on
Moritz [1972].  Collocation is a data processing method which
simultaneously perfoms regression, filtering and spatial prediction based
on given data.  Collocation methods generally assume that measurements
consist of a systematic component (trend surface) and a random component
(signal and noise):

B<w> = B<A> B<x> + B<s> + B<n>

with:

=over 4

=item B<A> Design matrix from constituent parts of mathematical model

=item B<x> Vector of unknown parameters (if any)

=item B<n> Vector of noise components in measurements

=item B<s> Vector of signal components in measurements

=item B<w> Vector of measurements

=back

The purpose of least-squares collocation is the prediction and the
separation of the random components B<s> and B<n> and furthermore,
the estimation of B<x> if required.

When the design matrix is known and removed from the measurements,
the task reduces to separate the signal and the noise components:

B<z> = B<w> - B<A> B<x> = B<s> + B<n>

The constitutents of the vector B<z> can be calculated from the
C<statistical quadratic form> which is the expression which
weighted least-squares adjustment would minimise:

sqf(B<v>) = B<v>^{T} Cov(B<v>,B<v>)^{-1} B<v>

The vectors B<v> contain the distances between the known samples
and the distances to the required point.  The covariance matrix
Cov(B<v>,B<v>) is the sum of the signal and noise covariance matrices.
The element I<ij> of the signal covariance matrix can be computed using
the Gaussian function:

I<sigs>(I<s>) = I<sigs>**2 * exp( - I<s>**2 / I<sh>**2 )

where I<s> is the separator factor representing the distance between
the points of interest, and I<sh> represents some constant control
parameter.  The signal variance <sigs> is set by the user.

If the observed measurements are error-free, then the noise
covariance matrix can be set to zero.  Otherwise a small error
is assigned to each available point:

I<sign>(I<s>) = I<sign>**2

The noise variance <sign> is set by the user.  By adjusting the
variances <sigs> and <sign> the user can control the separation
of the corresponding constituents.


=head2 EXPORT

Nothing.

=head1 METHODS

=head2 CONSTRUCTOR

=over 4

=item B<new>(I<key> => I<value>)

The constructor B<new> expects a list of key-value pairs and returns the
reference to the created instance.  The following keys are accepted.

=over 4

=item B<dim> => I<int>

The key B<dim> is required if you want to specify the dimension of the problem
only.  An instance is created without any samples.  Samples can be added later
using the method B<add> described below.

=item B<sample> => I<arrayref>

The key B<sample> can be used to submit existing data points of arbitrary
dimension.  If the submitted array reference contains numeric values only,
the problem is assumed to be 1D and the key B<value> is mandatory to submit
the corresponding values using a reference to an array of the same size.
If the submitted array reference contains arrays, then the problem is assumed
to be I<N>-dimensional according to the size of the largest array.
If the key B<value> is omitted, the values are taken from the last entry of
the available arrays and the dimension of the problem is reduced to I<N-1>.

=item B<value> => I<arrayref>

The key B<value> receives a reference to an array containing numeric values
that describe the properties to the corresponding samples given by the key
B<sample>.  The user has to take care that the sizes of the arrays referenced
by B<value> and B<sample> are the same.

=back

=cut

sub BUILD
{
    # Get caller reference
    my ($self) = @_;

    # Check arguments
    if (not defined $self->sample and not defined $self->dim)
    {
        # Neither option `sample` nor option `dim` defined
        croak "Method `new` requires either option `sample` or option `dim`";
    }
    elsif (not defined $self->sample)
    {
        # Option `dim` defined, but not `sample`
        croak "Method `new` requires positive dimension" if $self->dim < 1;

        # Clear points and values
        $self->{sample} = [];
        $self->{value} = [];

        # Return
        return;
    }

    # Flag to indicate that an array references is required for `sample`
    my $sample_need_arrayref = not defined $self->value;

    # Compare array sizes of `point` and `value`
    croak "Method `new` requires matching array sizes for `sample` and `value`"
        if $sample_need_arrayref and scalar @{$self->value} != $self->size;

    # Declare array for point dimensions
    my @pointdim;

    # Check array `sample` contains array references only
    for my $i (1 .. $self->size)
    {
        # Check array element and determine corresponding dimension
        my $elem = $self->{sample}[$i - 1];
        my $reftype = reftype($elem);
        my $dim = 0;
        if (not defined $reftype and not $sample_need_arrayref)
        {
            # Scalar found
            $dim = 1;
            $reftype = 'SCALAR';
        }
        elsif ($reftype eq 'ARRAY')
        {
            # Array found
            $dim = scalar @{$elem};
        }
        else
        {
            # Raise error if neither scalar nor array found
            croak "Method `new` requires array of arrays for key `sample`";
        }

        # Save detected dimension of current point
        $pointdim[$i - 1] = $dim;
    }

    # If option `value` is not given
    if ($sample_need_arrayref)
    {
        # Loop over elements of array `sample`
        for my $i (1 .. $self->size)
        {
            # Remove last element and shift it to array `value`
            my $value = pop @{$self->{sample}[$i - 1]};
            $self->{value}[$i - 1] = $value;
            $pointdim[$i - 1]--;
        }
    }

    # Check for unique dimensions
    my @uniqnum = sort { $a <=> $b } uniqnum @pointdim;

    # Require unique array elements when array references are required
    croak "Method `new` requires array sizes to be unique for key `sample`"
        if $sample_need_arrayref and scalar @uniqnum != 1;

    # Raise warning if maximum dimension of samples exceeds given dimension
    carp sprintf(
        "Method `new` requires dimension %d, but received `%d`",
            max(@uniqnum), $self->dim)
                if defined $self->dim and $self->dim < max(@uniqnum);

    # Adjust problem dimension (regardless if already defined)
    $self->dim(max(@uniqnum));

    # Add potentially mssing zeroes to each element of array `sample`
    for my $i (1 .. $self->size)
    {
        push @{$self->{sample}[$i - 1]}, 0
            while scalar @{$self->{sample}[$i - 1]} < $self->dim;
    }
}


=back

=head2 OBJECT METHODS

Here is a list of object methods available.  Object methods are applied to
the object in question, in contrast with class methods which are applied to
a class.

=over 4

=item B<size>

The method B<size> returns the number of available samples.

=cut

sub size
{
    my ($self) = @_;
    return scalar @{$self->sample};
}

=item B<sigs>(I<value>)

The method B<sigs> reads or sets the signal variance for the corresponding
object.  The default value is C<5>.


=item B<sign>(I<value>)

The method B<sign> reads or sets the noise variance for the corresponding
object.  Set the noise variance to zero, if the values at the given samples
are absolutely free of any measurement errors.  The default value is C<1>.

=item B<add>({ sample => I<arrayref>, value => I<arrayref> })

The method B<add> adds a single sample to the corresponding object.  The
method accepts a single argument, which can either be a hash reference or an
array reference.

=cut

sub add
{
    # Get caller reference and argument
    my ($self, $arg) = @_;
    my $reftype = reftype($arg);
    my @point;

    # Check reference type
    if (defined $reftype and $reftype eq 'HASH')
    {
        # Raise error, keys `sample` and `value` do not exist
        croak "Method `add` requires keys `sample` and `value`"
            if not exists $arg->{sample} and not exists $arg->{value};

        # Raise error, invalid type of `value`
        croak "Method `add` expects scalar for hash key `value`"
            if defined reftype($arg->{value});

        # Raise error, invalid type of `sample`
        croak "Method `add` expects array references for hash key `sample`"
            if $self->dim > 1 and
                (not defined reftype($arg->{sample}) or
                             reftype($arg->{sample}) ne 'ARRAY');

        # Check problem is 1D and `sample` is scalar
        if ($self->dim == 1 and not defined reftype($arg->{sample}))
        {
            # Add scalar to temporary array
            @point = ( $arg->{sample} );
        }
        else
        {
            # Dereference given array
            @point = @{$arg->{sample}};
        }

        # Add value to temporary array
        push @point, $arg->{value} if defined $arg->{value};
    }
    elsif (defined $reftype and $reftype eq 'ARRAY')
    {
        # Raise error, size or referenced array is invalid
        croak sprintf("Method `add` expects array of size %d, received: %d",
                $self->dim + 1, scalar @{$arg})
                    if $self->dim + 1 != scalar @{$arg};

        # Dereference given array
        @point = @{$arg};
    }
    else
    {
        # Raise error, argument has invalid type
        croak "Method `add` requires array or hash reference";
    }

    # Process temporary array
    if (all { defined $_ and m/^$RE{num}{real}$/i } @point)
    {
        # Remove last element
        my $value = pop @point;

        # Append data to keys `sample` and `value`
        push @{$self->{sample}}, [ @point ];
        push @{$self->{value}}, $value;
    }
    else
    {
        # Issue warning
        carp sprintf(
            "Method `add` ignores sample " .
            "containing non-numeric or undefined values: %s",
                ( join ',', map { defined $_ ? "`$_`" : 'undef' } @point ));
    }
}

sub _autorange
{
    # Get caller reference and dimension
    my ($self, $dim) = @_;

    # Check dimension if defined
    croak sprintf(
        "Method `_autorange` expects positive integer < %d", $self->dim)
            if defined $dim and $dim !~ /^$RE{num}{int}$/ and
                ($dim < 1  or $dim > $self->dim);

    # Init array reference on result
    my $result = [];

    # Check wheter distinct dimension given
    if (defined $dim)
    {
        # Return boundaries of corresponding dimension only
        $result = [ $self->_minval($dim), $self->_maxval($dim) ];
    }
    else
    {
        # Loop over all dimensions
        for my $i (1 .. $self->dim)
        {
            # Add boundaries to results array
            push @{$result}, [ $self->_minval($i), $self->_maxval($i) ];
        }
    }

    # Return array reference
    return $result;
}

sub _minval
{
    # Get caller reference and dimension
    my ($self, $dim) = @_;

    # Check dimension
    croak sprintf(
        "Method `_minval` expects positive integer < %d", $self->dim)
            if not defined $dim or $dim !~ /^$RE{num}{int}$/ or
                ($dim < 1  or $dim > $self->dim);

    # Return minimum value of corresponding dimension
    return min( map { $_->[ $dim - 1 ] } @{$self->{sample}} );
}

sub _maxval
{
    # Get caller reference and dimension
    my ($self, $dim) = @_;

    # Check dimension
    croak sprintf(
        "Method `_maxval` expects positive integer < %d", $self->dim)
            if not defined $dim or $dim !~ /^$RE{num}{int}$/ or
                ($dim < 1  or $dim > $self->dim);

    # Return maximum value of corresponding dimension
    return max( map { $_->[ $dim - 1 ] } @{$self->{sample}} );
}

sub _gaussian
{
    # Get caller reference and distance value
    my ($self, $value) = @_;

    # Return weighting factor
    return $self->sigs**2 * exp( -log(2) * $value**2 / $self->sh**2 );
}


=item B<prediction>

The method B<prediction> returns the normalized prediction fraction in the
range [0..1].

=cut

sub prediction
{
    my ($self) = @_;
    return $self->sigs**2 / ( $self->sigs**2 + $self->sign**2 );
}

=item B<filter>

The method B<filter> returns the normalized filtering fraction in the
range [0..1].

=cut

sub filter
{
    my ($self) = @_;
    return $self->sign**2 / ( $self->sigs**2 + $self->sign**2 );
}

=item B<interpolate>({ range => I<arrayref>, grid => I<arrayref> })

The method B<interpolate> applies the collocation algorithm and serves as
primary method of this class.  The method receives an optional argument which
can be of any type.  No argument is assumed as reference to an empty hash.
If the method is called in scalar context, the method will return a reference
to the array containing the results of the collocation.  If the method is
called in array context, the method will return the resulting array directly.


=over 4

=item B<hashref> (see example above)

If the submitted argument is a hash reference, then the keys B<range>
and B<grid> receive the boundaries and the coarseness of the domain of
interest.

=over 4

=item B<range> => I<arrayref>

The optional key B<range> receives an array reference either to scalars
indicating the lower and upper boundary of the requested range separated
by a colon (C<:>) or to arrays with 2 values, i.e. the lower and the
upper boudnary itselves.  If several array elements are undefined or the
array is not defined at all, the range is determined automatically based
on the minimum and the maximum value given by the samples.

=item B<grid> => I<int>|I<arrayref>

The optional key B<grid> may receive a positive integer which is applied
to all dimensions available or it may receive a array of integers which
indicate the requested grid size for each dimension.  If key is omitted,
the default grid of C<20> is applied.

=back

=item B<arrayref>

If the submitted value is an array reference, the value is assumed to be a
single sample with size corresponding to the dimension of the object.
The method will return the resulting value in scalar context.

=item B<scalar>

A scalar value is allowed for objects with dimension 1 only, otherwise an
error is raised.  The method will return the resulting value as a scalar.

=back

=cut

sub interpolate
{
    # Get caller reference and argument
    my ($self, $arg) = @_;

    # Check whether samples are available
    croak "Method `interpolate` cannot be applied, no data available"
        if $self->size == 0;

    # Declare local variables
    my ($result, $range, $grid, $reftype) = ( [], [], [] );

    # Check whether argument is submitted
    if (defined $arg)
    {
        # Get type of argument
        $reftype = reftype($arg);

        # Evaluate type
        if (defined $reftype and $reftype eq 'HASH')
        {
            # Type is a hash reference
            $range = $arg->{range} if exists $arg->{range};
            $grid  = $arg->{grid}  if exists $arg->{grid};

            # Error if value `grid` is invalid
            croak "Method `interpolate` expects scalar or array " .
                  "for argument `grid`"
                    if defined reftype($grid) and
                               reftype($grid) ne 'ARRAY';

            # Use default grid if undefined or empty array
            if (not defined $grid or
               (reftype($grid) eq 'ARRAY' and not @{$grid}))
            {
                # Raise warning and apply default
                carp "Method `interpolate` uses default grid";
                $grid = 20;
            }

            # Convert scalar to array
            $grid = [ $grid ] if not defined reftype($grid);

            # Raise error if array contains too many values
            croak sprintf(
                "Method `interpolate` expects %s value%s for `grid`, " .
                "received: %s",
                    $self->dim,
                    ( $self->dim > 1 ) ? 's' : '',
                    join(',', map { defined $_ ? "`$_`" : 'undef' } @{$grid}))
                        if scalar @{$grid} > $self->dim;

            # Add max grid to missing dimensions
            my $max = max(@{$grid});
            for my $i (1 .. $self->dim)
            {
                # Set to maximum if not defined
                $grid->[$i - 1] = $max unless defined $grid->[$i - 1];
            }

            # Raise error for non-positive integers
            croak sprintf(
                "Method `interpolate` received invalid element in `grid`: %s",
                    join(',', map { defined $_ ? "`$_`" : 'undef' } @{$grid}))
                        if not all { m/^$RE{num}{int}$/i and $_ >= 0 } @{$grid};
        }
        elsif (defined $reftype and $reftype eq 'ARRAY')
        {
            # Argument is an array
            croak sprintf(
                "Method `interpolate` expects array of size %d, received: %d",
                    $self->dim, scalar @{$arg})
                        if $self->dim != scalar @{$arg};

            # Adjust range to point
            for my $val (@{$arg})
            {
                push @{$range}, [ $val, $val ];
                push @{$grid},  [ 0 ];
            }
        }
        elsif (not defined $reftype and $self->dim == 1)
        {
            # Adjust range to point
            $range = [ [ $arg, $arg ] ];
            $grid  = [ [ 0 ] ];
        }
        else
        {
            # Raise error
            croak "Method `interpolate` expects array or hash reference";
        }
    }

    # Test range
    for my $i (1 .. $self->dim)
    {
        # Fetch element
        my $test = $range->[$i - 1];
        if (defined $test)
        {
            # Fetch element type
            if (defined reftype($test))
            {
                # Raise error if type is not an array
                croak "Method `interpolate` expects array reference"
                    unless reftype($test) eq 'ARRAY';

                # Raise error unless number of elements in array is not 2
                croak "Method `interpolate` expects array with 2 elements"
                    unless scalar @{$test} == 2;
            }
            else
            {
                # Element type is a scalar
                my ($lbound, $ubound) = split /\:/, $test;

                # Clear local array
                $test = [];

                # Check lower boundary
                if (not defined $lbound or $lbound !~ m/^$RE{num}{real}$/i)
                {
                    # Issue warning
                    carp sprintf(
                        "Method `interpolate` adjusts lower boundary of " .
                        "range for dimension %d", $i);

                    # Auto-determine lower boundary
                    $lbound = $self->_minval($i);
                }
                else
                {
                    # Assign detected value
                    $test->[0] = $lbound;
                }

                # Check upper boundary
                if (not defined $ubound or $ubound !~ m/^$RE{num}{real}$/i)
                {
                    # Issue warning
                    carp sprintf(
                        "Method `interpolate` adjusts upper boundary of " .
                        "range for dimension %d", $i);

                    # Auto-determine upper boundary
                    $ubound = $self->_maxval($i);
                }
                else
                {
                    # Assign detected value
                    $test->[1] = $ubound;
                }

                # Check which boundaries are user-defined
                if (not defined $test->[0] and not defined $test->[1])
                {
                    # Assign auto-detected values on both boundaries
                    $test = [ $lbound, $ubound ];
                }
                elsif (not defined $test->[0])
                {
                    # Assign auto-detected value on lower boundary
                    $test->[0] = min($lbound, $test->[1]);
                }
                elsif (not defined $test->[1])
                {
                    # Assign auto-detected value on upper boundary
                    $test->[1] = max($ubound, $test->[0]);
                }
            }

            # Check for correct range
            croak sprintf(
                "Method `interpolate` invoked with invalid boundaries " .
                "on dimension %d: %s",
                    $i,
                    join (',', map { defined $_ ? "`$_`" : 'undef' } @{$test}))
                        if not all { defined $_ and
                                     m/^$RE{num}{real}$/i } @{$test} or
                           $test->[0] > $test->[1];
        }
        else
        {
            # Auto-determine range
            $test = $self->_autorange($i);

            # Raise warning if auto-detection has been used
            carp sprintf(
                "Method `interpolate` misses range for dimension %d: " .
                "using %s", $i,
                    join(':', map { defined $_ ? "`$_`" : 'undef' } @{$test}));
        }

        # Apply local range
        $range->[$i - 1] = $test;
    }

    # Reset grid to zero where range is zero
    for my $i (1 .. $self->dim)
    {
        $grid->[$i - 1] = 0
            if isequal( $range->[$i - 1][1], $range->[$i - 1][0] );
    }

    # Submit grid to object
    $self->{grid} = $grid;

    # Number of available samples
    my $m = $self->size;

    # Check if reference length is defined
    my $shflag = defined $self->sh;

    # Calculate all spatial distances
    my $dist = Math::MatrixReal->new( $m, $m );

    # Loop over all points
    for my $k (1 .. $m)
    {
        # Point is zero distance from itself
        $dist->assign($k, $k, 0);

        # Get distances to all other points
        for my $l (($k + 1) .. $m)
        {
            # Calculate distance
            my $val = 0;
            for my $j (1 .. $self->dim)
            {
                $val += ( $self->{sample}[$k - 1][$j - 1]
                        - $self->{sample}[$l - 1][$j - 1] )**2;
            }
            $val = sqrt($val);

            # Check against minimum non-zero value
            $self->sh($val)
                if not $shflag and $val > 0
                               and (not defined $self->sh or $val < $self->sh);

            # Assign distance
            $dist->assign( $k, $l, $val );
            $dist->assign( $l, $k, $val );
        }
    }

    # Initialize signal and noise covariance matrices
    my $css = Math::MatrixReal->new( $m, $m );
    my $cnn = Math::MatrixReal->new( $m, $m );

    # Loop over all elements
    for my $k (1 .. $m)
    {
        # Assign diagonal of noise matrix
        $cnn->assign($k, $k, $self->sign**2 );

        #
        for my $l ($k .. $m)
        {
            # Fetch distance
            my $skl = $dist->element($k,$l);
            my $val = $self->_gaussian($skl);

            # Assign elements of signal matrix
            $css->assign($k,$l,$val);
            $css->assign($l,$k,$val) if $k != $l;
        }
    }

    # Calculate sum of signal and noise covariance matrix
    my $czz = $css + $cnn;

    # Calculate inverse of covariance matrix
    my $inv = $czz->inverse;

    # Determine trend, i.e. average value
    my $z0 = sum(@{$self->{value}}) / $m;

    # Calculate function vector
    my $zdash = Math::MatrixReal->new( $m, 1 );
    for my $i (1 .. $m)
    {
        $zdash->assign( $i, 1, $self->{value}[$i - 1] - $z0 );
    }

    # Set current point to start position
    my $point = [];
    for my $i (1 .. $self->dim) { $point->[$i - 1] = $range->[$i - 1][0]; }

    # Main cycle loop
    CYCLE: while (1)
    {
        # Initialize distance vector for current point
        my $ct = Math::MatrixReal->new( 1, $m );

        # Loop over elements of distance vector
        for my $k (1 .. $m)
        {
            # Calculate weight
            my $val = 0;
            for my $j (1 .. $self->dim)
            {
                $val += ( $self->{sample}[$k - 1][$j - 1] -
                                         $point->[$j - 1] )**2;
            }
            $val = sqrt($val);
            $val = $self->_gaussian($val);

            # Assign value to element
            $ct->assign( 1, $k, $val );
        }

        # Calculate result
        my $val = $z0 + $ct->multiply($inv)->multiply($zdash)->element(1,1);

        # Push coordinates of point and value to results' array
        push @{$result}, [ @{$point}, $val ];

        # Loop over dimensions to increment current point
        for my $j (1 .. $self->dim)
        {
            # Check if upper boundary of dimension is reached
            if (isequal($point->[$j - 1], $range->[$j - 1][1]))
            {
                # Final coordinate reached
                last CYCLE if $j == $self->dim;

                # Reset coordinate
                $point->[$j - 1] = $range->[$j - 1][0];
            }
            else
            {
                # Check if grid is needed
                if ($grid->[$j -1] > 0)
                {
                    # Increment current point
                    $point->[$j - 1] +=
                        ( $range->[$j - 1][1]
                        - $range->[$j - 1][0] ) / $grid->[$j - 1];
                }

                # Terminate loop
                last;
            }
        }
    }

    # Convert array to scalar if input argument is scalar and object is 1D
    $result = pop @{$result->[0]}
        if defined $arg and not defined $reftype and $self->dim == 1;

    # Return result according to context
    return (wantarray) ? @{$result} : $result;
}

=item B<newline>(I<index>)

The method B<newline> returns C<true> or C<false> whether the sample
referred to by the given I<index> closes a block, i.e. is a multiple of
the grid points along the first dimension.  This method is particularly
useful to insert a newline when the output data shall be forwarded
to B<gnuplot>.

=cut

sub newline
{
    my ($self, $i) = @_;
    return (defined $i and $self->dim > 1
                       and ($i + 0) % ($self->{grid}[0] + 1) == 0);
}

=back

=head1 SUBROUTINES

The following routines are provided.

=over 4

=item B<isequal>(I<value1>, I<value2>[, I<tol>])

This subroutine compares two floating-point numbers I<value1> and
I<value2> and checks whether they are equal within a reasonable
tolerance I<tol>.

=cut

sub isequal
{
    my ($x,$y,$tol) = @_;
    $tol = 1E-06 unless defined $tol;
    return abs( $x - $y ) / ( abs($x) + abs($y) + $tol ) < $tol;
}

=back

=cut

no Moose;
__PACKAGE__->meta->make_immutable;


=head1 REFERENCES

=over 4

=item Höpcke, W. (1980):

C<Fehlerlehre und Ausgleichsrechnung>.  De Gruyter, Pages 211 f.

=item Moritz, H. (1972):

C<Advanced least-squares methods>.  Report No. 175, Department of Geodetic Science, Ohio State University.  132 Pages.

=item Ruffhead, A. (1987):

C<An introduction to least-squares collocation>.  Survey Review, Volume 29, Number 224, Pages 85-94.

=back

=head1 ACKNOWLEDGEMENTS

My father Joachim Boljen for guiding me the way through the mists of statistical geodesy for filtering, prediction and modeling.

=head1 COPYRIGHT AND LICENSE

MIT License

Copyright (c) 2021 Matthias Boljen

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

=cut

1;
