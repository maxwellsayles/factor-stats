StaticLibrary(target='factor-stats',
              source=['factor-stats.cc'],
              CXXFLAGS=['-O3', '-Wall', '-Werror', '-std=c++11'],
              CPPPATH=['..'])
SConscript('squfof-parigp/SConscript')
SConscript('flint-factor/SConscript')
SConscript('pari-factor/SConscript')
SConscript('sspar-factor/SConscript')
