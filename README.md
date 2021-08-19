# Escher.jl

This package implements a single [Makie](https://makie.juliaplots.org/stable/index.html)
recipe called `escherplot` that plots maps of metabolic models resembling its [namesake
GUI](https://escher.github.io/#/). It's primary purpose is to facilitate the generation of
high quality metabolic maps from within Julia.

## Example
```
using Escher
using CairoMakie

escherplot("maps/core-map.json")
```

## Attributes
The `escherplot` recipe exposes a number of custom attributes that can be used to modify the
basic metabolic map figure.
