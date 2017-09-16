# gochimithequeutils

## Installation

```sh
$> go get github.com/tbellembois/gochimithequeutils/...
```

## Installation of python bindings

```sh
$> cd $GOPATH/src/github.com/tbellembois/gochimithequeutils
$> go generate

$> export GODEBUG=cgocheck=0
$> python
>>> import chimitheque as ch
>>> print(ch.__doc__)
Package chimitheque provides useful chemistry related utilities.

>>> ch.LinearToEmpiricalFormula("CH3")
'CH3'

```

Alternatively, if the `gopy` python module from [go-python/gopy](https://github.com/go-python/gopy)
is available from your `$PYTHONPATH`, one can do the following:

```py
>>> import gopy
>>> ch = gopy.load("github.com/tbellembois/gochimithequeutils/chimitheque", lang="cffi")
gopy> inferring package name...
gopy> loading 'github.com/tbellembois/gochimithequeutils/chimitheque'...
gopy> importing 'github.com/tbellembois/gochimithequeutils/chimitheque'

>>> print(ch.__doc__)
Package chimitheque provides useful chemistry related utilities.

>>> ch.LinearToEmpiricalFormula("CH3")
'CH3'
```
