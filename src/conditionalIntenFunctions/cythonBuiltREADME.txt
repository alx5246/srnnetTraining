A. Lonsberry
November 2016

README for folder cythonBuilt

DESCRIPTION

    In order to speed up the functionality of some of the functions and methods, we can use Cython. In particular I took
    one of the slowest functions that performs Stochasitic Gradient Descent (SGD) upon data to train the CIFs and
    then used Cython.

    I wrote up this code a while ago and have mostly forgotten how to use it. I accidently deleted the compiled code and
    had to at least figure out how to recompile it. Cython works by taking a .pyx file and compiling into a .c cython
    file, and then complied into a .pyd file that can be used directly in the code.

    In order to do this stuff I had to make a .pyx file, which stores the function I want to compile, and a
    blahblahblah_setup.py file to inform the computer how to compile.

    I then go to the cmd and cd to the ../cythonBuilt folder. From here I run
       > python cifUpdateCifCoeffsSDG_setup.py build_ext --compiler=msvc --inplace
    which will build all the stuff.

    Somehow under the hood, my msvc compiler is correctly selected, I am not sure how this happens. But it does work.

    Note if I need to recompile it seems to be okay to delete the folder /build and tis contents along with the .c and
    .pyd files as they will get rebuilt.

