### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 8a1e0b2a-f5d0-11ef-1234-0123456789ab
begin
	using Pkg
	Pkg.activate(mktempdir())
	Pkg.add(["Plots", "Distributions", "PlutoUI"])
	using Plots, Distributions, PlutoUI
end

# ╔═╡ 8a1e0b2b-f5d0-11ef-1234-0123456789ab
md"""
# 🎱 Powerball Expected Value Model

An interactive model for computing the expected value of a Powerball ticket as a function of the advertised jackpot, accounting for:

- **Lump-sum discount** (annuity → cash)
- **Federal + state taxes**
- **Ticket sales frenzy** (exponential growth with jackpot size)
- **Prize sharing** (Poisson model for multiple winners)
- **Non-jackpot prize EV**

Adjust the sliders below to explore different scenarios.
"""

# ╔═╡ 8a1e0b2c-f5d0-11ef-1234-0123456789ab
md"## Parameters"

# ╔═╡ 8a1e0b2d-f5d0-11ef-1234-0123456789ab
md"**Advertised Jackpot (millions USD):**"

# ╔═╡ 8a1e0b2e-f5d0-11ef-1234-0123456789ab
@bind J Slider(100:10:3000, default=800, show_value=true)

# ╔═╡ 8a1e0b2f-f5d0-11ef-1234-0123456789ab
md"**Frenzy coefficient** (higher = more ticket sales at high jackpots):"

# ╔═╡ 8a1e0b30-f5d0-11ef-1234-0123456789ab
@bind c Slider(0.0005:0.0001:0.0070, default=0.0030, show_value=true)

# ╔═╡ 8a1e0b31-f5d0-11ef-1234-0123456789ab
md"**Lump-sum factor** (fraction of advertised jackpot paid as cash):"

# ╔═╡ 8a1e0b32-f5d0-11ef-1234-0123456789ab
@bind lumpSumFactor Slider(0.40:0.01:0.65, default=0.50, show_value=true)

# ╔═╡ 8a1e0b33-f5d0-11ef-1234-0123456789ab
md"**Tax retention factor** (fraction kept after federal + state taxes):"

# ╔═╡ 8a1e0b34-f5d0-11ef-1234-0123456789ab
@bind taxFactor Slider(0.50:0.01:0.80, default=0.63, show_value=true)

# ╔═╡ 8a1e0b35-f5d0-11ef-1234-0123456789ab
md"**Non-jackpot prize EV** (expected value from smaller prizes, per ticket):"

# ╔═╡ 8a1e0b36-f5d0-11ef-1234-0123456789ab
@bind otherPrizesEV Slider(0.20:0.01:0.60, default=0.40, show_value=true)

# ╔═╡ 8a1e0b37-f5d0-11ef-1234-0123456789ab
md"**Baseline ticket sales** `a` (millions, at low jackpots):"

# ╔═╡ 8a1e0b38-f5d0-11ef-1234-0123456789ab
@bind a Slider(1.0:0.1:20.0, default=5.0, show_value=true)

# ╔═╡ 8a1e0b39-f5d0-11ef-1234-0123456789ab
md"**Sales scale factor** `b` (millions):"

# ╔═╡ 8a1e0b3a-f5d0-11ef-1234-0123456789ab
@bind b Slider(1.0:0.5:50.0, default=10.0, show_value=true)

# ╔═╡ 8a1e0b3b-f5d0-11ef-1234-0123456789ab
md"## Model"

# ╔═╡ 8a1e0b3c-f5d0-11ef-1234-0123456789ab
md"""
**Combinatorics:** Powerball has ``\binom{69}{5} \times 26 = 292{,}201{,}338`` possible tickets.

**Ticket sales model:**
```math
n(J) = a + b \, e^{c \, J}
```
where ``J`` is the advertised jackpot in millions, ``a`` is baseline sales, ``b`` is a scale factor, and ``c`` is the frenzy coefficient.

**Expected winners (Poisson parameter):**
```math
\lambda = \frac{n \times 10^6}{292{,}201{,}338}
```

**Prize sharing factor:** If ``K \sim \text{Poisson}(\lambda)`` other tickets also win, your expected share is:
```math
S(\lambda) = \sum_{k=0}^{\infty} \frac{1}{k+1} \, P(K = k) = E\!\left[\frac{1}{K+1}\right]
```

**Expected value per ticket:**
```math
\text{EV} = \underbrace{f_{\text{lump}} \cdot f_{\text{tax}}}_{\text{after-tax lump sum}} \cdot \underbrace{\frac{1}{292{,}201{,}338}}_{\text{win probability}} \cdot J \times 10^6 \cdot S(\lambda) \;+\; \text{EV}_{\text{other prizes}}
```
"""

# ╔═╡ 8a1e0b3d-f5d0-11ef-1234-0123456789ab
begin
	# Constants
	const COMBOS = 292_201_338
	const P_WIN = 1.0 / COMBOS
	const TICKET_PRICE = 2.0
end;

# ╔═╡ 8a1e0b3e-f5d0-11ef-1234-0123456789ab
"""
    share_factor(λ; max_k=150)

Compute E[1/(K+1)] where K ~ Poisson(λ), truncating the sum at max_k.
"""
function share_factor(λ; max_k=150)
	kmax = min(max_k, ceil(Int, λ) + 30)
	return sum((1.0 / (k + 1)) * pdf(Poisson(λ), k) for k in 0:kmax)
end

# ╔═╡ 8a1e0b3f-f5d0-11ef-1234-0123456789ab
"""
    compute_ev(J, a, b, c, lumpSumFactor, taxFactor, otherPrizesEV)

Compute the after-tax, lump-sum-adjusted expected value of a Powerball ticket.
"""
function compute_ev(J, a, b, c, lumpSumFactor, taxFactor, otherPrizesEV)
	ntickets = a + b * exp(c * J)          # millions of tickets sold
	λ = ntickets * 1e6 / COMBOS            # expected number of winners
	sf = share_factor(λ)                    # E[1/(K+1)]
	jackpot_ev = P_WIN * (J * 1e6) * sf    # raw jackpot contribution
	return lumpSumFactor * taxFactor * jackpot_ev + otherPrizesEV
end

# ╔═╡ 8a1e0b40-f5d0-11ef-1234-0123456789ab
begin
	ntickets = a + b * exp(c * J)
	λ = ntickets * 1e6 / COMBOS
	sf = share_factor(λ)
	jackpot_ev = P_WIN * (J * 1e6) * sf
	full_ev = lumpSumFactor * taxFactor * jackpot_ev + otherPrizesEV
end;

# ╔═╡ 8a1e0b41-f5d0-11ef-1234-0123456789ab
md"""
## Results

| Quantity | Value |
|:---------|------:|
| Tickets sold | **$(round(ntickets, digits=1)) million** |
| Expected winners (λ) | **$(round(λ, digits=4))** |
| Share factor S(λ) | **$(round(sf, digits=5))** |
| Raw jackpot EV | **\$$(round(jackpot_ev, digits=3))** |
| After-tax lump-sum jackpot EV | **\$$(round(lumpSumFactor * taxFactor * jackpot_ev, digits=3))** |
| Non-jackpot prize EV | **\$$(round(otherPrizesEV, digits=2))** |
| **Total EV per \$2 ticket** | **\$$(round(full_ev, digits=3))** |
| **Net EV (profit per ticket)** | **\$$(round(full_ev - TICKET_PRICE, digits=3))** |

$(full_ev > TICKET_PRICE ? "✅ **Positive expected value!**" : "❌ **Negative expected value** — the house wins on average.")
"""

# ╔═╡ 8a1e0b42-f5d0-11ef-1234-0123456789ab
md"## EV vs. Jackpot Size"

# ╔═╡ 8a1e0b43-f5d0-11ef-1234-0123456789ab
let
	jj_range = 100:25:3000
	ev_values = [compute_ev(jj, a, b, c, lumpSumFactor, taxFactor, otherPrizesEV)
	             for jj in jj_range]

	plot(jj_range, ev_values,
		label="Adjusted EV",
		xlabel="Advertised Jackpot (millions USD)",
		ylabel="EV per \$2 ticket (USD)",
		title="Powerball Expected Value vs. Jackpot Size",
		linewidth=3,
		color=:steelblue,
		legend=:topright,
		size=(800, 500),
		dpi=150,
		margin=5Plots.mm)

	hline!([TICKET_PRICE], label="Breakeven (\$2)", linestyle=:dash, color=:red, linewidth=2)

	# Mark current jackpot
	vline!([J], label="Current J = \$$(J)M", linestyle=:dot, color=:gray, linewidth=1.5)
	scatter!([J], [full_ev], label="", markersize=8, color=:orange)
end

# ╔═╡ 8a1e0b44-f5d0-11ef-1234-0123456789ab
md"## Sensitivity: Frenzy Coefficient"

# ╔═╡ 8a1e0b45-f5d0-11ef-1234-0123456789ab
let
	jj_range = 100:25:3000
	cs = [0.001, 0.002, 0.003, 0.005, 0.007]

	p = plot(xlabel="Advertised Jackpot (millions USD)",
		ylabel="EV per \$2 ticket (USD)",
		title="EV Sensitivity to Frenzy Coefficient c",
		legend=:topright,
		size=(800, 500),
		dpi=150,
		margin=5Plots.mm)

	for ci in cs
		ev_vals = [compute_ev(jj, a, b, ci, lumpSumFactor, taxFactor, otherPrizesEV)
		           for jj in jj_range]
		plot!(p, jj_range, ev_vals, label="c = $ci", linewidth=2)
	end

	hline!([TICKET_PRICE], label="Breakeven", linestyle=:dash, color=:red, linewidth=2)
	p
end

# ╔═╡ 8a1e0b46-f5d0-11ef-1234-0123456789ab
md"""
## Key Insights

1. **The frenzy kills the value.** As jackpots grow, ticket sales explode exponentially, driving up λ and destroying the share factor. The EV curve *peaks and then declines* at very high jackpots.

2. **Breakeven is rare.** Under realistic parameters (c ≈ 0.003, lump-sum ≈ 50%, tax retention ≈ 63%), EV almost never exceeds \$2.

3. **The only path to positive EV** is weak frenzy (low c), which would require a jackpot that somehow doesn't generate excitement — a contradiction in practice.

4. **Non-jackpot prizes** contribute a fixed ~\$0.40 per ticket regardless of jackpot size, providing a floor on total EV.
"""

# ╔═╡ Cell order:
# ╟─8a1e0b2a-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b2b-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b2c-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b2d-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b2e-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b2f-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b30-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b31-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b32-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b33-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b34-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b35-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b36-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b37-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b38-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b39-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b3a-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b3b-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b3c-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b3d-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b3e-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b3f-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b40-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b41-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b42-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b43-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b44-f5d0-11ef-1234-0123456789ab
# ╠═8a1e0b45-f5d0-11ef-1234-0123456789ab
# ╟─8a1e0b46-f5d0-11ef-1234-0123456789ab
