import sys
import os


# this variable should contain the names of all libraries we build
our_libs = [
    'lse_vtk',
    'mls_vtk',
    'mss_vtk'
    ]


from info import build_dir, win32_build_type

prefix = ''
if os.name == 'posix':
    sys.path.append(os.path.join(build_dir, 'bin'))
    prefix = 'lib'
elif os.name == 'nt':
    # we supose on windows visual studio is used to compile
    sys.path.append(os.path.join(build_dir, 'bin', win32_build_type))
    

for lib_name in our_libs:
    lib = __import__(prefix + lib_name + 'Python',
                     globals(),
                     locals())
    for k, v in lib.__dict__.items():
        if k[0] != '_':
            globals()[k] = v
