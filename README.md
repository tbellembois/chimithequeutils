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
