from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.command.build_ext import build_ext as du_build_ext


class build_ext(du_build_ext):
    def run(self):
        from Cython.Build.Dependencies import cythonize
        self.distribution.ext_modules[:] = cythonize(
            self.distribution.ext_modules,
            language_level=3)
        du_build_ext.run(self)


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

extensions = [
    Extension(
        "tencoloredpoints.tencoloredpoints",
        sources=["tencoloredpoints/tencoloredpoints.pyx"])
    ]

setup(
    name='tencoloredpoints',
    version='0.1.0',
    description='Show that each 10 points in the plane' +
            ' there exists a colored Tverberg partition',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/kliem/TenColoredPoints',
    author='Jonathan Kliem',
    author_email='jonathan.kliem@fu-berlin.de',
    license='GPLv3',
    packages=find_packages(),
    ext_modules=extensions,
    zip_safe=False,
    python_requires='>=3.6',
    package_dir={'tencoloredpoints': 'tencoloredpoints'},
    install_requires=["Cython", "memory_allocator"],
    package_data={"tencoloredpoints": ["*.pxd"]},
    cmdclass={'build_ext': build_ext},
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Mathematics']
    )
