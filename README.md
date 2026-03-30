# Powerball Expected Value Model

An interactive [Pluto.jl](https://plutojl.org/) notebook that computes the expected value of a Powerball lottery ticket as a function of the advertised jackpot size.

## Features

- **Interactive sliders** for all model parameters (jackpot, frenzy coefficient, lump-sum factor, tax rate, etc.)
- **Live-updating plots** showing EV vs. jackpot size
- **Sensitivity analysis** across different frenzy coefficients
- **Full mathematical derivation** embedded in the notebook

## The Model

The key insight: as jackpots grow, ticket sales explode exponentially, increasing the probability of splitting the prize. This "frenzy effect" means EV peaks at moderate jackpots and *declines* at very high ones.

**Ticket sales:** `n(J) = a + b * exp(c * J)`

**Expected winners:** `λ = n * 10⁶ / 292,201,338`

**Share factor:** `S(λ) = E[1/(K+1)]` where `K ~ Poisson(λ)`

**EV per ticket:** `f_lump * f_tax * P(win) * J * S(λ) + EV_other`

## Running Locally

1. Install [Julia](https://julialang.org/downloads/) (v1.10+)
2. Launch Pluto:
   ```julia
   using Pkg; Pkg.add("Pluto"); using Pluto; Pluto.run()
   ```
3. Open `powerball_ev.jl` in the Pluto interface

## Running in the Browser (No Install)

Click the badge below to launch an interactive session via Binder:

[![Open in Binder](https://mybinder.org/badge_logo.svg)](https://binder.plutojl.org/open?url=https://raw.githubusercontent.com/pyrex41/powerball-ev/main/powerball_ev.jl)

## License

MIT
