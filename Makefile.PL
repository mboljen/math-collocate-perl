use 5.006;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME          => 'Math::Collocate',
    AUTHOR        => q{Matthias Boljen <mboljen@online>},
    VERSION_FROM  => 'lib/Math/Collocate.pm',
    ABSTRACT_FROM => 'lib/Math/Collocate.pm',
    LICENSE       => 'MIT',
    PL_FILES      => {},
    EXE_FILES     => [],
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => 0,
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    PREREQ_PM => {
        'File::Basename'            => 0,
        'File::Temp'                => 0,
        'Getopt::Long'              => 0,
        'List::Util'                => 0,
        'Math::MatrixReal'          => 0,
        'Pod::Usage'                => 0,
        'Moose'                     => 0,
        'Scalar::Util'              => 0,
        'Spreadsheet::Read::Simple' => 0,
        'Time::HiRes'               => 0,
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'Math-Collocate-*' },
);
