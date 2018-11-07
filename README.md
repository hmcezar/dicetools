# DICEtools
Repository of selected Python 3 scripts used to aid data analysis and input generation of Monte Carlo and Configurational Bias Monte Carlo simulations performed with Dice.


## Dependencies
To use the scripts a few Python 3 dependencies are needed, for example
- [NumPy](http://www.numpy.org/)
- [OpenBabel](http://openbabel.org/)
- [matplotlib](https://matplotlib.org/)

The number of dependencies may vary depending on which analysis script one is using, those are just the most common used in almost all the scripts.
You can install these dependencies in any way you want, however, we encourage the use of the [Anaconda](https://www.anaconda.com/download/) Python distribution.
To install the libraries with Anaconda, do the following:
```
conda install numpy
conda install -c openbabel openbabel
conda install matplotlib
```

If you want to use the graphical interface of DiceWin to perform some analysis, you will also need some stuff like:
- [SciPy](https://www.scipy.org/) 
- [SIP](https://pypi.org/project/SIP/)
- [PyQt5](https://riverbankcomputing.com/software/pyqt/intro)
that can also be easily installed with conda with
```
conda install scipy
conda install sip
conda install pyqt=5
```


## Short description of tools
All the scripts can be run with the `-h` option to show a brief description of what the script does and the mandatory and optional parameters.


## Authorship
Most of the scripts here were written by Henrique Musseli Cezar, with the exception of DiceWin which was written by Thiago de Souza Duarte.
These tools were written with the important contribution of Prof. Kaline Coutinho, who supervised the work and gave suggestions to the improvement of the tools.


## Acknowledgments
We thank the Brazilian funding agencies [CNPq](http://www.cnpq.br/) and [CAPES](http://www.capes.gov.br/) for the fellowships and the approved research projects.
Most of this work was done under a CNPq PhD fellowship for Henrique Musseli Cezar (grant number 140489/2015-0) and undergraduate research fellowships for Thiago de Souza Duarte.