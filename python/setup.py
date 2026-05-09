from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [
    Extension ("pyatnept",
              ['pyatnept.pyx', 'PythonStream.cpp'],
              language = 'c++',
              extra_compile_args = ["-std=c++11"],
              libraries = ['atnept'],
              include_dirs = ['../planet', '../lib', '../tinyxml2'],
              library_dirs = ['..'],
              extra_link_args = ['-fopenmp'],
              )]

setup(
    name = 'pyatnept',
    ext_modules = cythonize(extensions, language_level = "2"),
    )
