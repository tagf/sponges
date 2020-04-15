# Sponges and totally incompressible porous medium

[Derivation of formulas using Python Sympy](/sympy_eqs_potential.ipynb)

[Our numerical simulation of 1D sponge using Python Numpy](/sponges_numeric.ipynb)

To reproduce the results, you might want to install Anaconda https://docs.anaconda.com/anaconda/install/
(automatically includes Python3, Spyder 4, and IPython Notebook support with JupyterNotebook / JupyterLab) on your computer.

Alternatively you can try JupyterLab online: https://jupyter.org/try or https://colab.research.google.com/notebooks/intro.ipynb

Here are direct links to open the notebooks above using google colab:

Formulas: https://colab.research.google.com/github/tagf/sponges/blob/master/sympy_eqs_potential.ipynb

!Unfortunately, google colab currently does not support rendering of LaTeX formulas in computation cells (I'm continuing to check if there's a way to do that) and it's default sympy version is is not up to date. To resolve the latter issue please uncomment (remove #  from) this line:

"# !pip install sympy==1.5.1"

Simulation:  https://colab.research.google.com/github/tagf/sponges/blob/master/sponges_numeric.ipynb


If you use JupyterLab online, click File, then create a new notebook,
paste the code from one of our computation .ipynb into a code cell of the notebook and run it.
