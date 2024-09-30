project('HeRoam', 'c')

# Source Files
srcs = files(
     'src/main.c',
     'src/vectors.c',
     'src/simulation.c'
)

# Build Type options
debug_flags = ['-g', '-pg', '-Wall', '-Wextra', '-fopenmp']
optimized_flags = ['-O3', '-march=native', '-Wall', '-Wextra', '-fopenmp', '-fopt-info-vec-optimized', '-fmax-errors=5']

# Flags based on buildtype
if get_option('buildtype') == 'debug'
   cflags = debug_flags
else
    cflags = optimized_flags
endif

# Linker flags
debug_ldflags = ['-pg', '-lm', '-fopenmp']
optimized_ldflags = ['-lm', '-fopenmp']

# Flags based on buildtype
if get_option('buildtype') == 'debug'
   ldflags = debug_ldflags
else
    ldflags = optimized_ldflags
endif


# Build
executable('HeRoam',
            srcs,
            c_args: cflags,
            link_args: ldflags)
