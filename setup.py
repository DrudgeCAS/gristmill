"""Setup script for gristmill."""

import os.path
import sys
from setuptools import setup, find_packages, Extension


PROJ_ROOT = os.path.dirname(os.path.abspath(__file__))
INCLUDE_DIRS = [
    os.path.join(PROJ_ROOT, 'deps', i, 'include')
    for i in ['cpypp', 'fbitset', 'libparenth']
]

# Platform-specific compiler flags
if sys.platform == "win32":
    # MSVC compiler flags
    COMPILE_FLAGS = ['/std:c++20']
else:
    # GCC/Clang compiler flags  
    COMPILE_FLAGS = ['-std=c++20']
    
    # Additional flags for macOS to avoid header conflicts
    if sys.platform == "darwin":
        # Use libc++ standard library explicitly on macOS
        COMPILE_FLAGS.extend([
            '-stdlib=libc++',
            '-mmacosx-version-min=10.9'
        ])

parenth = Extension(
    'gristmill._parenth',
    ['gristmill/_parenth.cpp'],
    include_dirs=INCLUDE_DIRS,
    extra_compile_args=COMPILE_FLAGS
)

setup(
    packages=find_packages(),
    ext_modules=[parenth],
    package_data={'gristmill': ['templates/*']},
)
