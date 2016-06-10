# voronoi

A testbed for a Voronoi diagram approach to monodromy and homology.

**Idea:** compute the Voronoi diagram on the branch point structure of a plane
algebraic curve. Using these data we can efficiently compute the monodromy group
of the curve and, from there, a basis for the first homology group. This
approach allows one to maximally re-use analtyic continuation information.

## Testing

To run tests:

```
$ sage runtests.py
```

To run tests with pdb (for interactive error debugging):

``` abap
$ sage -sh
(sage-sh) $ nosetests --pdb
```

## Authors

* Nils Bruin
* Sasha Zotine (@szotine)
* Chris Swierczewski (@cswiercz)
