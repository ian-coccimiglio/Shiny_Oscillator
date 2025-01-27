---
title: Simulating Biological Oscillations in R
author: Ian Coccimiglio
date: Sun 26 Jan 2025 17:41:51 PST

---

# Creating a Shiny app for interacting with biological oscillations

## Usage
Biological circuits are a fascinating way of understanding how something with complex dynamics can arise from relatively simple interactions. However, even in these simple situations, many parameters contribute to how the system evolves. As such, I simulated multiple systems of ordinary differential equations (ODEs) in R, and use cOde to compile them into C code for heightened performance (h/t to Thomas Petzoldt for this concept!). I then integrated these into a Shiny app to visualize 4 specific of creating a biological oscillator.

1) A damped biological oscillator, which has only two components, and is unable to oscillate permanently on its own. However, adding noise is one way to achieve indefinite (if unpredictable) oscillations.
2) A brusselator is a noise-free method of inducing biological oscillations with only 2 participating species, but requires additional reactions to achieve oscillations
3) Repressilation is a mechanism of inducing oscillations by having 3 candidate species all acting in repression, but requires specific parameterization of reaction rates to achieve stable oscillations.
4) Repressilation (with modifiable cooperativity) is a way of seeing the effect of cooperative binding on oscillatory mechanics, in which increasing cooperativity induces stavke oscillatory behaviour.

## Requirements:
deSolve
cOde
bslib
shiny

Tested on R version 4.40

## Acknowledgments
Thomas Petzoldt, !(for his work on translating ODEs written in R into compilable C code)[https://tpetzoldt.github.io/deSolve-shiny/deSolve-shiny.html#the\_brusselator\_in\_c]
Michael Elowitz, !(for his python code version of the repressilation system)[https://biocircuits.github.io/chapters/09\_repressilator.html]
Uri Alon, !(for his great book on biological circuits and systems biology)[https://www.amazon.ca/Introduction-Systems-Biology-Principles-Biological/dp/1439837171]
