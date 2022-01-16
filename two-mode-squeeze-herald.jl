### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 637715ed-ceb6-4992-b665-026088559278
begin
	using QuantumOptics, Plots, PlutoUI
	const F₀ = 6.5e3/2/π # 6.5 GHz omega
	const k_per_h = 20.84e3 # MHz per K
	const n_per_T = k_per_h/F₀ # photons per K
	const tsteps = 50
	
	# Raniwala's design
	const gᵐˢ = 12.
	const gᵒᵐ₀ = 0.268 # MHz
	const T₁ᵃ = 3e4 / 190e6 # 1/MHz
	const T₂ˢ = 1000. # Fake large number for an unused part of the system

	function prep(; N=3,
	                gᵒᵐ=0.268,
	                gᵐˢ = 12.,
	                T₁ᵃ = 3e4 / 190e6,
	                T₂ˢ = 1000.,
	                γ = 0.201449,
			        nₜ = 0.00325
	                )
	
	    ℬₚ = FockBasis(N)
	    ℬₛ = SpinBasis(1//2)
	    ℬ = ℬₚ⊗ℬₚ⊗ℬₛ
	
	    â = b̂ = destroy(ℬₚ)
	    n̂ᵃ = n̂ᵇ = number(ℬₚ)
	    σ̂₋ = sigmam(ℬₛ)
	    n̂ˢ = projector(spinup(ℬₛ))
	
	    Â   = embed(ℬ, 1, â)
	    B̂   = embed(ℬ, 2, b̂)
	    N̂ᵃ  = embed(ℬ, 1, n̂ᵃ)
	    N̂ᵇ  = embed(ℬ, 2, n̂ᵇ)
	    Σ̂₋  = embed(ℬ, 3, σ̂₋)
	    N̂ˢ  = embed(ℬ, 3, n̂ˢ)
	
	    Ĥᵒᵐ = gᵒᵐ * Â * B̂
	    Ĥᵒᵐ += Ĥᵒᵐ'
	
	    Ĥᵐˢ = gᵐˢ * B̂ * Σ̂₋'
	    Ĥᵐˢ += Ĥᵐˢ'
	
	    Ĉˢ = 1/√T₂ˢ * N̂ˢ
	    Ĉᵃ = 1/√T₁ᵃ * Â

	    Ĉᵇ = √(γ/2 * (nₜ+1)) * B̂
	    Ĉᵇ₁ = √(γ/2 * nₜ)    * B̂'
	
	    Ĥᵢₘ = -1im/2 * (Ĉᵃ'*Ĉᵃ + Ĉᵇ'*Ĉᵇ + Ĉᵇ₁'*Ĉᵇ₁ + Ĉˢ'*Ĉˢ)
	
	    vac  = fockstate(ℬₚ, 0) ⊗ fockstate(ℬₚ, 0) ⊗ spindown(ℬₛ)
	    vacꜛ = fockstate(ℬₚ, 0) ⊗ fockstate(ℬₚ, 0) ⊗ spinup(ℬₛ)
		excitationpair = fockstate(ℬₚ, 1) ⊗ fockstate(ℬₚ, 1) ⊗ spindown(ℬₛ)
	
	    Pvac  = projector(vac)
	    Pvacꜛ = projector(vacꜛ)
	    return (;(Base.@locals)...)
	end
	
	function apply_jump(ψ,c)
	    ψₙ = c*ψ
		ψₙ /= norm(ψₙ)
		ψₙ
	end
	
	pdf(ψ,c) = norm(c*ψ)/norm(ψ)
	
	md"Some more helper code in this cell."
end

# ╔═╡ b576bacc-c769-4b0d-b033-ad1ad0cd4181
md"""
# Heralded Single Phonon Generation in an Optomechanical System

_Supplementary to "Spin-optomechanical quantum interface enabled by an ultrasmall mechanical and optical mode volume cavity" by H. Raniwala, S. Krastanov, M. Eichenfield, and D. R. Englund, 2022_

Can be found at:

- Repository [`github.com/Krastanov/optomechspin_heralded_entanglement`](https://github.com/Krastanov/optomechspin_heralded_entanglement)
- Interactive at [`pluto.krastanov.org/two-mode-squeeze-herald.html`](https://pluto.krastanov.org/two-mode-squeeze-herald.html)
- Archived together with the publication

For any system in which two-mode squeezing is used for heralded generation of single photons.

$\hat{H}=g \hat{a}^\dagger \hat{b}^\dagger + H.c.
	- \frac{i}{2}\sum \hat{c}_l^\dagger\hat{c}_l$

with

- photon decay $\hat{c}_a=\frac{1}{\sqrt{T_a}}\hat{a}$
- phonon decay $\hat{c}_b=\sqrt{\gamma\frac{n_\textrm{th}+1}{2}}\hat{b}$
- phonon heating $\hat{c}_{b1}=\sqrt{\gamma\frac{n_\textrm{th}}{2}}\hat{b}^\dagger$

resulting into

$\hat{H}=g \hat{a}^\dagger \hat{b}^\dagger + H.c.
    -\frac{i}{2}\left(
		\frac{1}{T_a}n̂_a + γn_\textrm{th}\left(n̂_b+\frac{1}{2}\hat{Id}\right) 
	\right)$

We will be starting with $|00\rangle$. Most of the time no event will be happening, however an $a$ jump might happen leading to the desired heralded state. Leading sources of infidelity will be events $b$ or $b_1$ happening before or after $a$.

In the basis $\left\{| 00\rangle,| 11\rangle,| 01\rangle,| 10\rangle\right\}$ we get

$-i\hat{H}=
\begin{pmatrix}
0   & -ig                                        &   &   \\
-ig & -\frac{1}{2}(\frac{1}{T_a}+γ(n_\textrm{th}+½)) &   &   \\
    &                                            & -\frac{γ(n_\textrm{th}+½)}{2} & 0 \\
    &                                            & 0 & -\frac{1}{2T_a} 
\end{pmatrix}-\frac{\gamma n_\textrm{th}\hat{Id}}{4}$

**We neglect double excitations by choosing this basis.** We could have used a smaller basis $\left\{| 00\rangle,| 11\rangle\right\}$, but the larger one is instructive for when we look at the heralded jump.

In the smaller basis we have (with more convenient constants as defined below)

$-i\hat{H}=
\begin{pmatrix}
-\lambda_h   & -ig  \\
-ig & -2\lambda_l-\lambda_h
\end{pmatrix}$

The eigen solutions are

$|\psi_-\rangle = e^{-\lambda t} e^{-iGt} \frac{1}{\sqrt2}
\begin{pmatrix}
\frac{g}{G-i\lambda_l} \\
1
\end{pmatrix}$

and

$|\psi_+\rangle = e^{-\lambda t} e^{iGt} \frac{1}{\sqrt2}
\begin{pmatrix}
-\frac{g}{G+i\lambda_l} \\
1
\end{pmatrix}$

where

- total amplitude shrinkage $\lambda = \lambda_\textrm{l}+\lambda_\textrm{h} = \frac{1}{4T_a}+\frac{1}{4}γ(2n_\textrm{th}+½)$
- loss param $\lambda_\textrm{l} = \frac{1}{4}(\frac{1}{T_a}+γ(n_\textrm{th}+½))$
- heat param $\lambda_\textrm{h} = \frac{\gamma n_\textrm{th}}{4}$
- corrected interaction frequency $G = \sqrt{g^2-\lambda_l^2}$
- notice that $\left|\frac{g}{G+i\lambda_l}\right|=\left|\frac{g}{G-i\lambda_l}\right|=1$ hence the $\frac{1}{\sqrt2}$ normalization factor
"""

# ╔═╡ 74b628b2-82d5-471b-9f40-18b3b547668b
md"""
## Evolution starting at $|00\rangle$

For initial state $| 00\rangle$ this implies

$|\psi_{00}\rangle = \frac{g}{\sqrt2G}\left(|\psi_-\rangle - |\psi_+\rangle\right)
= e^{-\lambda t}
\begin{pmatrix}
\cos(Gt) + \frac{\lambda_l}{G}\sin(Gt) \\
-i \frac{g}{G} \sin(Gt)
\end{pmatrix}$

## Evolution starting at $|11\rangle$

$|\psi_{11}\rangle =
\frac{1}{\sqrt2G}\left(\left(G-i\lambda_l\right)|\psi_-\rangle + \left(G+i\lambda_l\right)|\psi_+\rangle\right)
= e^{-\lambda t}
\begin{pmatrix}
-i\frac{g}{G}\sin(Gt) \\
\cos(Gt) - \frac{\lambda_l}{G}\sin(Gt)
\end{pmatrix}$
"""

# ╔═╡ 89a7d7af-34a5-434d-b409-b40a070b5a3c
md"""
## What if $\lambda_l > g$?

Then $G$ becomes imaginary, but taking into account that $i\sinh(x)=\sin(ix)$ and $\cosh(x)=\cos(ix)$ we get

$|\psi_{00}\rangle = e^{-\lambda t}
\begin{pmatrix}
\cosh(|G|t) + \frac{\lambda_l}{|G|}\sinh(|G|t) \\
-i \frac{g}{|G|} \sinh(|G|t)
\end{pmatrix}
\underset{t\gg 1/\lambda_l,\ g\ll\lambda_l}{\to}
\frac{1}{2}e^{-\left(\frac{g^2}{2\lambda_l} +\lambda_h\right)t}
\begin{pmatrix}
1 + \frac{\lambda_l}{|G|} \\
-i \frac{g}{|G|}
\end{pmatrix}$

$|\psi_{11}\rangle = e^{-\lambda t}
\begin{pmatrix}
-i\frac{g}{|G|}\sinh(|G|t) \\
\cosh(|G|t) - \frac{\lambda_l}{|G|}\sinh(|G|t)
\end{pmatrix}
\underset{t\gg 1/\lambda_l,\ g\ll\lambda_l}{\to}
\frac{1}{2}e^{-\left(\frac{g^2}{2\lambda_l} +\lambda_h\right)t}
\begin{pmatrix}
-i\frac{g}{|G|} \\
1 - \frac{\lambda_l}{|G|}
\end{pmatrix}$

**Notice that if $g\ll\lambda_l$, then the effect of the $|G|-\lambda_l$ in the exponential is quadratically suppressed and a quasi-steady state is reached very rapidly and then decays at rate $\frac{g^2}{\lambda}$.**
"""

# ╔═╡ a7d208cc-5c47-4fb3-8920-472e224bdc53
md"""

### Comparison to exact calculation

Below you can compare our approximate analytic solution to a complete numerical simulation and see they match as long as we do not keep state $|11\rangle$ populated (in order to avoid double excitations).

`g :` $(@bind gtoy Slider(0.1:0.1:5; default=1, show_value=true))

`λₗ:` $(@bind λₗtoy Slider(0.1:0.1:5; default=2, show_value=true))

`λₕ:` $(@bind λₕtoy Slider(0.0:0.1:5; default=0, show_value=true))

Log axis $(@bind logaxistoy CheckBox())
"""

# ╔═╡ aba7b5fd-9905-4915-a691-34cee7fbbb7b
md"""## State after jump at time $\tau$

Various jumps can occur and we will need to be able to track them if we are to calculate infidelities of the heralded state (as opposed to calculating just the heralding probability).

- after $\hat{c}_a$: $| \psi_a \rangle = | 01\rangle$
- after $\hat{c}_b$: $| \psi_b \rangle = | 10\rangle$
- after $\hat{c}_{b1}$: $| \psi_{b1} \rangle =$ superposition of $| 01\rangle$ and $| 12 \rangle$

$\hat{H}=g \hat{a}^\dagger \hat{b}^\dagger + H.c.
    -\frac{i}{2}\left(
		\frac{1}{T_a}n̂_a + γn_\textrm{th}\left(n̂_b+\frac{1}{2}\hat{Id}\right) 
	\right)$

In the basis $\left\{| 00\rangle,| 11\rangle,| 01\rangle,| 10\rangle,| 12\rangle,| 21\rangle\right\}$ we get

$-i\hat{H}=
\begin{pmatrix}
0   & -ig                                        &   &  & & \\
-ig & -\frac{1}{2}(\frac{1}{T_a}+γ(n_\textrm{th}+½)) &   &  & & \\
    &                                            & -\frac{γ(n_\textrm{th}+½)}{2} & 0& -\sqrt2ig & \\
    &                                            & 0 & -\frac{1}{2T_a} & 0&-\sqrt2ig \\
&&-\sqrt2ig& 0 & -\frac{1}{2T_a} -γ(n_\textrm{th}+½) & 0 \\
&&&-\sqrt2ig& 0 & -\frac{1}{T_a} -\frac{γ(n_\textrm{th}+½)}{2}
\end{pmatrix}-\frac{\gamma n_\textrm{th}\hat{Id}}{4}$

Conveniently, we can again restrict ourselves to two-dimensional subspaces, either $\left\{| 01\rangle, | 12\rangle\right\}$, or $\left\{| 10\rangle, | 21\rangle\right\}$, in both of which we have a Hamiltonian with the same shape as previously solved, but different values for the parameters:

$-i\hat{H}=
\begin{pmatrix}
-\tilde{\lambda}_h   & -i\tilde{g}  \\
-i\tilde{g} & -2\tilde{\lambda}_l-\tilde{\lambda}_h
\end{pmatrix}$

**We are again neglecting double excitations by choosing this basis.**

### Dynamics of $|01\rangle$

In $\left\{| 01\rangle, | 12\rangle\right\}$ we have

- Interaction rate $\tilde{g} = \sqrt2 g$
- Excited state decay param $\tilde{\lambda}_l = \lambda_l = \frac{1}{4}(\frac{1}{T_a}+γ(n_\textrm{th}+½))$ 
- Shared state decay param $\tilde{\lambda}_{h01} = \gamma\left(\frac{3}{4}n_\textrm{th}+\frac{1}{4}\right)$

### Dynamics of $|10\rangle$

In $\left\{| 10\rangle, | 21\rangle\right\}$ we have

- Interaction rate $\tilde{g} = \sqrt2 g$
- Excited state decay param: $\tilde{\lambda}_l = \lambda_l = \frac{1}{4}(\frac{1}{T_a}+γ(n_\textrm{th}+½))$ 
- Shared state decay param $\tilde{\lambda}_{h10} = \frac{1}{2T_a}+\frac{1}{4}\gamma n_\textrm{th}$

"""

# ╔═╡ 2fdaf940-174b-411f-905f-7d3dd101dc2f
begin
	G(g,λₗ) = √Complex(g^2-λₗ^2)
	λₗ(Tₐ, γ, nₜₕ) = 0.25*(1/Tₐ + γ*(nₜₕ+0.5))
    λₕ(Tₐ, γ, nₜₕ) = 0.25*γ*nₜₕ
    λₕ₀₁(Tₐ, γ, nₜₕ) = γ*(0.75*nₜₕ+0.25)
    λₕ₁₀(Tₐ, γ, nₜₕ) = 0.5/Tₐ + 0.25*γ*nₜₕ
	"""The time evolution operator, i.e. columns {|ψ₀₀⟩,|ψ₁₁⟩}"""
	function U(g,λₗ,λₕ,t)
		Γ = G(g,λₗ)
		λ = λₗ+λₕ
		c = cos(Γ * t) 
		s = sin(Γ * t)
		exp(-λ*t)*[c+λₗ/Γ*s   -im*g/Γ*s;
		           -im*g/Γ*s  c-λₗ/Γ*s]
	end
	"""The time evolution operator in the g≪λₗ limit."""
	function Uₗᵢₘ(g,λₗ,λₕ,t)
		Γ = abs(G(g,λₗ))
		0.5 * exp((-g^2/2/λₗ + λₕ)*t) * [ 1+λₗ/Γ -im*g/Γ;
										 -im*g/Γ 1-λₗ/Γ]
	end
	md"Implementing the aforementioned formulas"
end

# ╔═╡ e6b3962c-1a45-47ac-ad12-0e429c2bcccc
let
	ts0 = (0:0.025:4)
	ts = ts0./gtoy # time span
	# Theory plot
	Us = cat(U.(gtoy,λₗtoy,λₕtoy,ts)...,dims=3)
	Usapprox = cat(Uₗᵢₘ.(gtoy,λₗtoy,λₕtoy,ts)...,dims=3)
	p00 = plot(ts0, abs.(Usapprox[1,1,:]), label=nothing, c=1,
		title="Amplitudes vs time for different initial states")
	plot!(ts0, abs.(Usapprox[1,2,:]), label=nothing, c=2)
	plot!(ts0, abs.(Us[1,1,:]), label=nothing, c=1, ls=:dot)
	plot!(ts0, abs.(Us[1,2,:]), label=nothing, c=2, ls=:dot)
	p11 = plot(ts0, abs.(Usapprox[2,1,:]), label="|00⟩ (g≪λₗ limit)", c=1,
		xlabel="time (1/g)")
	plot!(ts0, abs.(Usapprox[2,2,:]), label="|11⟩ (g≪λₗ limit)", c=2)
	plot!(ts0, abs.(Us[2,1,:]), label="|00⟩ (two-level approx)", c=1, ls=:dot)
	plot!(ts0, abs.(Us[2,2,:]), label="|11⟩ (two-level approx)", c=2, ls=:dot)
	# Numerics plot
	nₜ = 1.
	γ = 4*λₕtoy/nₜ
	T₁ᵃ = 1/4/(λₗtoy+γ*(nₜ+0.5))
	e = prep(;gᵒᵐ=gtoy,gᵐˢ=0,T₁ᵃ=T₁ᵃ,T₂ˢ=1e6,nₜ=nₜ,γ=γ)
	# TODO find a way to use the nice DifferentialEquations.jl interface
	_, ψ₀₀ = timeevolution.schroedinger(ts, e.vac, e.Ĥᵒᵐ+e.Ĥᵢₘ)
	ψ₀₀00 = (ψ -> e.vac'*ψ).(ψ₀₀)
	ψ₀₀11 = (ψ -> e.excitationpair'*ψ).(ψ₀₀)
	plot!(p00,ts0,abs.(ψ₀₀00),ls=:dash,c=1, label=nothing)
	plot!(p00,ts0,abs.(ψ₀₀11),ls=:dash,c=2, label=nothing)
	_, ψ₁₁ = timeevolution.schroedinger(ts, e.excitationpair, e.Ĥᵒᵐ+e.Ĥᵢₘ)
	ψ₁₁00 = (ψ -> e.vac'*ψ).(ψ₁₁)
	ψ₁₁11 = (ψ -> e.excitationpair'*ψ).(ψ₁₁)
	plot!(p11,ts0,abs.(ψ₁₁00),ls=:dash,c=1, label="|00⟩ (full numerics)")
	plot!(p11,ts0,abs.(ψ₁₁11),ls=:dash,c=2, label="|11⟩ (full numerics)")

	if logaxistoy
		plot(p00,p11,layout=(2,1),yaxis=:log10,ylims=(1e-3,1))
	else
		plot(p00,p11,layout=(2,1))
	end
end

# ╔═╡ 8f5f0424-1e03-473e-8cd4-32934fadc581
md"""

## Probability of event

The probability density for an event $c$ is

$\textrm{pdf}_c = \frac{\textrm{d}P_c}{\textrm{d}t}
= \langle\psi\mid\hat{c}^\dagger\hat{c}\mid\psi\rangle
/ \langle\psi\mid\psi\rangle$

and the probability of any event happening at all is

$\int_0^t\textrm{pdf}_* = \int_0^t\sum_i \textrm{pdf}_{c_i}=1-\langle\psi(t)\mid\psi(t)\rangle.$ 

Here are the branches of the multiverse that we need to track (events that "click" as heralded, but only the first one is good):

- **just $a$ (the good event)**: In the limit of large $t$ and small $g$ its PDF is $\textrm{pdf}_a(t) = \frac{1}{T_a}\left(\frac{g}{2\lambda_l}\right)^2$.
- **an $a$ at $t=\tau$ followed by any other event**: Creates an $|01\rangle$ and evolves in $\{|01\rangle,|12\rangle\}$. In the same limit its conditional probability after the $a$ jump is $P_{a*}(t,\tau)=1-e^{-\left(\frac{\tilde{g}^2}{\lambda_l}+2\tilde{\lambda}_{h01}\right)(t-\tau)}$
- **a $b$ at $t=\tau$ followed by $a$**: Creates an $|10\rangle$ and evolves in $\{|10\rangle,|21\rangle\}$. In the same limit its probability density for the $b$ jump is $\textrm{pdf}_b(t) = \gamma\frac{n_\textrm{th}+1}{2}\left(\frac{g}{2\lambda_l}\right)^2$ and the probability denstity for $a$ after $b$ is $\textrm{pdf}_{ba}(t_{>\tau})=\frac{1}{T_a}$
- **a $b_1$ followed by $a$**: Creates a superposition evolving in $\{|01\rangle,|12\rangle\}$. In the same limit its probability density for the $b_1$ jump is $\textrm{pdf}_{b1}(t) = \gamma\frac{n_\textrm{th}}{2}$ and the probability denstity for $a$ after $b_1$ is $\textrm{pdf}_{b1a}(t_{>\tau})=\frac{1}{T_a}\left(\frac{\tilde{g}}{2\lambda_l}\right)^2$

Thus, for pumps of duration $T$, as long as the henceforth calculated probabilities are small, we end up with:

- Probability of $a$ happening first $P_\textrm{first a} = \int_0^T \textrm{d}t\ \textrm{pdf}_a(t)  = \frac{T}{T_a}\left(\frac{g}{2\lambda_l}\right)^2$
- Probability of $a$ happening first, followed by a bad event $P_\textrm{a then bad}=\int_0^T\textrm{d}\tau\ \textrm{pdf}_a(\tau) P_{a*}(T,\tau) =  \frac{T+C^{-1}(e^{-CT}-1)}{T_a}\left(\frac{g}{2\lambda_l}\right)^2$, with $C=\frac{\tilde{g}^2}{\lambda_l}+2\tilde{\lambda}_{h01}$
- Probability of $b$ followed by $a$ $P_{ba} = \int_0^T \textrm{d}\tau\ \textrm{pdf}_b(\tau) \int_\tau^T\textrm{d}t\ \textrm{pdf}_{ba}(t)=\frac{T^2}{2}\frac{1}{T_a}\gamma\frac{n_\textrm{th}+1}{2}\left(\frac{g}{2\lambda_l}\right)^2$
- Similarly for $b_1$ followed by $a$ $P_{b1a} =\frac{T^2}{2}\frac{1}{T_a}\gamma\frac{n_\textrm{th}}{2}\left(\frac{\tilde{g}}{2\lambda_l}\right)^2$
"""

# ╔═╡ 5ea35653-e67d-4f33-9bbd-2f59d1c731a7
begin
	pdf_a(t,g,λₗ,Tₐ) = 1/Tₐ * (g/2/λₗ)^2
	P_a✽(t,τ,g,λₗ,λₕ₀₁) = 1 - exp(-2* (g^2/λₗ + λₕ₀₁) * (t-τ))
	pdf_b(t,g,λₗ,γ,nₜₕ) = γ * (nₜₕ+1)/2 * (g/2/λₗ)^2
	pdf_ba(t,Tₐ) = 1/Tₐ
	pdf_b1(t,γ,nₜₕ) = γ*nₜₕ/2
	pdf_b1a(t,g,λₗ,Tₐ) = pdf_a(t,√2*g,λₗ,Tₐ)
	md"Implementing the aforementioned formulas"
end

# ╔═╡ 05bcbadf-b764-4f0a-8115-8938860794ea
md"""

### Comparison to exact calculation

Comparing the probability densities derived above to numerical simulations. This time we will use values based on the design by Hamza Raniwala.

Pump photons
`log₁₀(n) : `$(@bind log₁₀n Slider(0:0.5:3; default=2, show_value=true))

Temperature (in K)
`log₁₀(T) : `$(@bind log₁₀T Slider(-2:0.1:2; show_value=true))

Mechanical Quality Factor
`log₁₀(Qₘ): `$(@bind log₁₀Qₘ Slider(3:0.5:9; default=6, show_value=true))

Pump Duration (in units of optical lifetime)
`τ = t/T₁ᵃ: `$(@bind τ Slider(1:100; default=20, show_value=true))

Jump location (if any)
`jump indx: `$(@bind i_jmp Slider(2:tsteps-2; default=tsteps÷2))
"""

# ╔═╡ a03f56e0-bf39-4645-a57f-0955903cf7a2
begin
	n = 10. ^ log₁₀n
	gᵒᵐ = √n * gᵒᵐ₀ # MHz
	Temp = 10. ^ log₁₀T
	Qₘ = 10. ^ log₁₀Qₘ
	nₜ = n_per_T*Temp
	γ = π*F₀/Qₘ

	md"""
	from Raniwala's design:
	$g_0$ = $(round(gᵒᵐ₀,digits=3)) MHz and T_a = $(round(T₁ᵃ*1000,digits=3)) ns
	(i.e. 1/Tₐ= $(round(1/T₁ᵃ/1000,digits=2)) GHz)
	
	derived:
	
	| n | g | T | Qₘ | nₜₕ |
	|:--|:--|:--|:---|:----|
	| $(round(n)) | $(round(gᵒᵐ,digits=1)) MHz | $(round(Temp,digits=3)) K | $(round(Qₘ,digits=-3)) | $(round(nₜ,digits=2))
	"""
end

# ╔═╡ 69ad33a1-7497-4439-bdfa-8a8e30bd5068
let
	Λₗ = λₗ(T₁ᵃ,γ,nₜ)
	Λₕ₀₁ = λₕ₀₁(T₁ᵃ,γ,nₜ)
	
	e = prep(;gᵒᵐ,gᵐˢ,T₁ᵃ,T₂ˢ,nₜ,γ)

	# Before Jump
	times = range(0, τ*T₁ᵃ; length=tsteps)
	Δt = times[2]
	_, states = timeevolution.schroedinger(times, e.vac, e.Ĥᵒᵐ+e.Ĥᵢₘ)
	# TODO: tracking the spin is not actually necessary here
	
	# After Jump
	t_jmp = times[i_jmp]
	times_ajmp = range(0, τ*T₁ᵃ*(tsteps-i_jmp)/tsteps; length=tsteps-i_jmp)
	state_ajmp = e.Ĉᵃ*states[i_jmp]
	state_ajmp /= norm(state_ajmp)
	_, states_ajmp = timeevolution.schroedinger(times_ajmp, state_ajmp, e.Ĥᵒᵐ+e.Ĥᵢₘ)
	state_bjmp = e.Ĉᵇ*states[i_jmp]
	state_bjmp /= norm(state_bjmp)
	_, states_bjmp = timeevolution.schroedinger(times_ajmp, state_bjmp, e.Ĥᵒᵐ+e.Ĥᵢₘ)
	state_b1jmp = e.Ĉᵇ₁*states[i_jmp]
	state_b1jmp /= norm(state_b1jmp)
	_, states_b1jmp = timeevolution.schroedinger(times_ajmp, state_b1jmp, e.Ĥᵒᵐ+e.Ĥᵢₘ)
	
	make_pdf(Op) = ψ -> real(expect(Op'*Op , ψ))/norm(ψ)
	make_pop(Op) = ψ -> real(expect(Op , ψ))/norm(ψ)
	
	C̄ᵃ  = make_pdf(e.Ĉᵃ ).(states)
	C̄ᵇ  = make_pdf(e.Ĉᵇ ).(states)
	C̄ᵇ₁ = make_pdf(e.Ĉᵇ₁).(states)
	C̄ˢ  = make_pdf(e.Ĉˢ ).(states)
	plot( times, C̄ᵃ , label="a",  c=1, lw=2)
	plot!(times, C̄ᵇ , label="b",  c=2)
	plot!(times, C̄ᵇ₁, label="b₁", c=3)
	#plot!(times, C̄ˢ , label="s",  c=4)
	C̄ᵃ_ajmp  = make_pdf(e.Ĉᵃ ).(states_ajmp)
	C̄ᵇ_ajmp  = make_pdf(e.Ĉᵇ ).(states_ajmp)
	C̄ᵇ₁_ajmp = make_pdf(e.Ĉᵇ₁).(states_ajmp)
	C̄ˢ_ajmp  = make_pdf(e.Ĉˢ ).(states_ajmp)
	plot!(times_ajmp.+t_jmp,
		C̄ᵃ_ajmp+C̄ᵇ_ajmp+C̄ᵇ₁_ajmp+C̄ˢ_ajmp,c=4,label="any jump after a")
	# Too large, but conditioned on too small
	#C̄ᵃ_bjmp  = make_pdf(e.Ĉᵃ ).(states_bjmp)
	#plot!(times_ajmp.+t_jmp,
	#	C̄ᵃ_bjmp,c=5,label="a jump after b")
	C̄ᵃ_b1jmp  = make_pdf(e.Ĉᵃ ).(states_b1jmp)
	plot!(times_ajmp.+t_jmp,
		C̄ᵃ_b1jmp,c=5,label="a jump after b1")
	plot!(times, pdf_a.(times, gᵒᵐ,Λₗ,T₁ᵃ), ls=:dash, label="when g≪λₗ", c=1)
	plot!(times, pdf_b.(times, gᵒᵐ,Λₗ,γ,nₜ), ls=:dash, label=nothing, c=2)
	plot!(times, pdf_b1.(times,γ,nₜ), ls=:dash, label=nothing, c=3)
	pdf_a✽ = P_a✽.(τ*T₁ᵃ,t_jmp,gᵒᵐ,Λₗ,Λₕ₀₁) ./ (τ*T₁ᵃ-t_jmp)
	pdf_a✽ = pdf_a✽ * ones(length(times_ajmp))
	plot!(times_ajmp.+t_jmp, pdf_a✽,c=4, label=nothing, ls=:dash)
	# Too large, but conditioned on too small
	#plot!(times_ajmp.+t_jmp, pdf_ba.(times_ajmp,T₁ᵃ),c=5, label=nothing, ls=:dash)
	plot!(times_ajmp.+t_jmp, pdf_b1a.(times_ajmp,gᵒᵐ,Λₗ,T₁ᵃ),c=5, label=nothing, ls=:dash)
	p_pdf_before = plot!(title="Probability density for event",
		xlabel="time", legend=:topleft)

	#plot(times_ajmp.+t_jmp , C̄ᵃ_ajmp , label=nothing, c=1, lw=2)
	#plot!(times_ajmp.+t_jmp , C̄ᵇ_ajmp , label=nothing, c=2)
	#plot!(times_ajmp.+t_jmp , C̄ᵇ₁_ajmp, label=nothing, c=3)
	#plot!(times_ajmp.+t_jmp , C̄ˢ_ajmp , label=nothing, ls=:dash, c=4)
	#p_pdf_after = plot!(link=:all, xlabel="time")
	#plot(p_pdf_before, p_pdf_after, layout=(2,1))
    	
	N̄ᵃ  = make_pop(e.N̂ᵃ ).(states)
	N̄ᵇ  = make_pop(e.N̂ᵇ ).(states)
	N̄ˢ  = make_pop(e.N̂ˢ ).(states)
	plot( times, N̄ᵃ , label="a",  c=1)
	plot!(times, N̄ᵇ , label="b",  c=2)
	#plot!(times, N̄ˢ , label="s",  c=4)
	p_pop = plot!(title="Population", xlabel="time")

	N̄ᵃ_ajmp  = make_pop(e.N̂ᵃ ).(states_ajmp)
	N̄ᵇ_ajmp  = make_pop(e.N̂ᵇ ).(states_ajmp)
	N̄ˢ_ajmp  = make_pop(e.N̂ˢ ).(states_ajmp)
	plot( times_ajmp.+t_jmp, N̄ᵃ_ajmp , label="a",  c=1, ls=:dash, lw=2)
	plot!(times_ajmp.+t_jmp, N̄ᵇ_ajmp , label="b",  c=2, ls=:dash)
	#plot!(times_ajmp.+t_jmp, N̄ˢ_ajmp , label="s",  c=4, ls=:dash)
	p_pop_ajmp = plot!(title="Population after jump", xlabel="time", xlim=(0,times[end]))

	plot( times, cumsum(C̄ᵃ)*Δt , label="a", lw=2)
	plot!(times, cumsum(C̄ᵇ)*Δt , label="b")
	plot!(times, cumsum(C̄ᵇ₁)*Δt, label="b₁")
	#plot!(times, cumsum(C̄ˢ)*Δt , label="s")
	plot!(times, cumsum(C̄ᵃ+C̄ᵇ+C̄ᵇ₁+C̄ˢ)*Δt , label="sum", ls=:dash, lc=:gray)
	plot!(times, 1 .- norm.(states).^2, label="1-norm²", lc=:black)
	p_norm = plot!(title="Cumulative chance of event", xlabel="time")
	
	#plot(p_pdf_before,p_pop, p_pop_ajmp, p_norm, layout=(4,1), size=(600,700))
	plot(p_pdf_before, p_norm, layout=(2,1), size=(600,700), xlabel="time (μs)")
end

# ╔═╡ 13430563-5293-4484-8b37-fcab0e442bd4
md"""
## Heralding Probability and Fidelity

Finally, we get probability of a heralding click (**we assume the detector can not distinguish two $a$ photons in rapid succession from a single $a$ photon**)

$P_\textrm{herald} = P_\textrm{first a} + P_{ba} + P_{b1a} =
\frac{T}{T_a}\left(\frac{g}{2\lambda_l}\right)^2 \left(1+T\frac{\gamma}{2}\frac{3n_\textrm{th}+1}{2}\right)$

The state being heralded is

$\rho = p_\textrm{good\ traj}\ \rho_\textrm{good\ traj} + p_\textrm{bad}\rho_\textrm{bad}$

where $p_\textrm{good\ traj} = \frac{P_\textrm{first a} - P_\textrm{a then bad}}{P_\textrm{herald}} = \frac{1-e^{-CT}}{CT}\frac{1}{1+T\frac{\gamma}{2}\frac{3n_\textrm{th}+1}{2}}$ and $p_\textrm{bad} = 1-p_\textrm{good\ traj}$.

Thus the final fidelity is

$F =\langle 1|\rho|1\rangle= p_\textrm{good\ traj} \langle 1|\rho_\textrm{good\ traj}|1\rangle
= \frac{1-e^{-CT}}{CT}\frac{1}{1+T\frac{\gamma}{2}\frac{3n_\textrm{th}+1}{2}} \langle 1|\rho_\textrm{good\ traj}|1\rangle.$

**We have pesimistically assumed $\langle 1|\rho_\textrm{bad}|1\rangle=0$. We also have an implicit trace with respect to the photonic mode somewhere in here.**
"""

# ╔═╡ 87b86d27-520c-4947-91cb-ce321f245d26
md"""
## The drift of the heralded state

We need to find $\rho_\textrm{good\ traj}$, the state averaged over all trajectories that have a single good heralding event and no bad events. We just need to average over the probability density for such an event.

$\rho_\textrm{good\ traj} = \left.\int_0^T\textrm{d}t\ \textrm{pdf}_a(t)\rho_\textrm{after\ jump}(T-t) \middle/ \int_0^T\textrm{d}t\ \textrm{pdf}_a(t)\right.$

where $\rho_\textrm{after\ jump}(T-t)$ is the state $|01\rangle$ evolved for time $T-t$ and then partial-traced with respect to the photonic mode:

$\rho_\textrm{after\ jump}(T-t)=|\alpha|^2 |1\rangle\langle 1| +|\beta|^2|2\rangle\langle 2|$

$\alpha = \frac{1}{2}e^{-\left(\frac{\tilde g^2}{2\lambda_l} +\tilde \lambda_{h01}\right)t}\left(1 + \frac{\lambda_l}{|\tilde G|}\right)$

$\beta = \frac{1}{2}e^{-\left(\frac{\tilde g^2}{2\lambda_l} +\tilde \lambda_{h01}\right)t}\left(-i \frac{\tilde g}{|\tilde G|}\right)$

Leading to 

$\langle 1|\rho_\textrm{good\ traj}|1\rangle = \frac{1-e^{-CT}}{CT}$
"""

# ╔═╡ 9c0a28b6-fb44-46cc-85d1-4c1daa4e098a
md"""
## Complete expression for the fidelity and success probability

$P =
\frac{T}{T_a}\left(\frac{g}{2\lambda_l}\right)^2 \left(1+T\frac{\gamma}{2}\frac{3n_\textrm{th}+1}{2}\right)$

$F = \left(\frac{1-e^{-CT}}{CT}\right)^2\frac{1}{1+T\frac{\gamma}{2}\frac{3n_\textrm{th}+1}{2}}$

$C=\frac{\tilde{g}^2}{\lambda_l}+2\tilde{\lambda}_{h01}$

$\tilde{g} = \sqrt2 g$
$\lambda_l = \frac{1}{4}(\frac{1}{T_a}+γ(n_\textrm{th}+½))$ 
$\tilde{\lambda}_{h01} = \gamma\left(\frac{3}{4}n_\textrm{th}+\frac{1}{4}\right)$

And after plugging everything in:

$P =
\frac{T}{T_a}\left(\frac{g}{\frac{1}{2}(\frac{1}{T_a}+γ(n_\textrm{th}+½))}\right)^2
\left(1+\frac{\gamma T}{2}\frac{3n_\textrm{th}+1}{2}\right)$

$F = \left(\frac{1-e^{-CT}}{CT}\right)^2\frac{1}{1+\frac{\gamma T}{2}\frac{3n_\textrm{th}+1}{2}}$

$C=
\frac{2g^2}{\frac{1}{4}(\frac{1}{T_a}+γ(n_\textrm{th}+½))}
+\gamma\frac{3n_\textrm{th}+1}{2}$

## Design rules of thumb

The probability $P$ is monotonicly increasing and the fidelity $F$ is monotonicly decreasing with $T$, so for high fidelity we would just pick the minimal permitted $T$ at which the steady state under consideration is possible, i.e., $T\gg \left(T_a^{-1} + \left(\frac{1}{γ(n_\textrm{th}+½)}\right)^{-1}\right)^{-1}$. This is easy to achive while keeping $F$ high only if $T_a \ll T \ll \frac{1}{γ(n_\textrm{th}+½)}$. Therefore the regime of interest is:

$T_a \ll T \ll \frac{1}{γ(n_\textrm{th}+½)}$
$T_a \ll \frac{1}{g}$

On the other hand, we want to maximize the probability of success $P$, or rate of success $P/T$ under these same constraints, which requires a high $g$.

Thus, in this protocol, first and foremost is the optimization of the mechanical Q factor (or its thermal population).
"""

# ╔═╡ 72621af3-af5b-412d-86cb-a6376f983a5b
md"""
In the aforementioned parameter regime we get:

$C\approx
8g^2T_a
+\gamma\frac{3n_\textrm{th}+1}{2}$

and 

$P \approx 4 g^2 T_a T$

$1 - F \approx 1 - \frac{1-CT}{1+\frac{\gamma T}{2}\frac{3n_\textrm{th}+1}{2}} \approx 8g^2T_aT+\frac{3}{4}\gamma T(3n_\textrm{th}+1)$


"""

# ╔═╡ 5b6c3ce3-adfe-414d-b74f-65f4084f49ac
begin
	P(T,g,Tₐ,γ,nₜₕ) = T/Tₐ * (g/2/λₗ(Tₐ,γ,nₜₕ))^2 * (1+γ*T/4*(3nₜₕ+1))
    F(T,g,Tₐ,γ,nₜₕ) = (1-exp(-C(g,Tₐ,γ,nₜₕ)*T))^2 / (C(g,Tₐ,γ,nₜₕ)*T)^2 / (1+γ*T/4*(3nₜₕ+1))
	C(g,Tₐ,γ,nₜₕ) = 8*g^2 / (1/Tₐ+γ*(nₜₕ+1/2))  +  γ/2*(3nₜₕ+1)
	Tminapprox(g,Tₐ,γ,nₜₕ) = 1/(1/Tₐ+γ*(nₜₕ+1/2))
	Papprox(T,g,Tₐ,γ,nₜₕ) = 4*g^2*T*Tₐ
	Fapprox(T,g,Tₐ,γ,nₜₕ) = 1 - 8*g^2*T*Tₐ - 3/4*γ*T*(3nₜₕ+1)
	md"Implementing the aforementioned formulas"
end

# ╔═╡ 35667cfb-b837-4197-b7d5-f59334109969
md"""

Pump photons
`log₁₀(n) : `$(@bind log₁₀n2 Slider(0:0.5:3; default=2, show_value=true))

Temperature (in K)
`log₁₀(T) : `$(@bind log₁₀T2 Slider(-2:0.1:2; show_value=true))

Mechanical Quality Factor
`log₁₀(Qₘ): `$(@bind log₁₀Qₘ2 Slider(3:0.5:9; default=6, show_value=true))
"""

# ╔═╡ 4ff08761-cc25-4e04-9849-492332c39fd8
begin
	n2 = 10. ^ log₁₀n2
	#gᵒᵐ₀ = 0.268 # MHz defined above
	gᵒᵐ2 = √n2 * gᵒᵐ₀ # MHz
	#T₁ᵃ = 3e4 / 190e6 # 1/MHz defined above
	Temp2 = 10. ^ log₁₀T2
	Qₘ2 = 10. ^ log₁₀Qₘ2
	nₜ2 = n_per_T*Temp2
	γ2 = π*F₀/Qₘ2

	md"""
	from Raniwala's design:
	$g_0$ = $(round(gᵒᵐ₀,digits=3)) MHz and T_a = $(round(T₁ᵃ*1000,digits=3)) ns
	(i.e. 1/Tₐ= $(round(1/T₁ᵃ/1000,digits=2)) GHz)
	
	derived:
	
	| n | g | T | Qₘ | nₜₕ |
	|:--|:--|:--|:---|:----|
	| $(round(n2)) | $(round(gᵒᵐ2,digits=1)) MHz | $(round(Temp2,digits=3)) K | $(round(Qₘ2,digits=-3)) | $(round(nₜ2,digits=2))
	"""
end

# ╔═╡ 78f9535c-9952-4b80-b681-73f9496371b1
let
	ts = exp.(range(log(0.1),log(100),length=100)).*T₁ᵃ
	minT = Tminapprox(gᵒᵐ2,T₁ᵃ,γ2,nₜ2)
	ts = ts[ts .> minT]
	ps = P.(ts,gᵒᵐ2,T₁ᵃ,γ2,nₜ2)
	infs = 1 .- F.(ts,gᵒᵐ2,T₁ᵃ,γ2,nₜ2)
	plot(infs,ps,line_z=log10.(ts),lw=4,
		xlabel="1-F", ylabel="P", colorbar_title=" \nlog₁₀(T/μs)",
		right_margin=3Plots.PlotMeasures.mm,
		label=nothing,
		xscale=:log10,yscale=:log10,
		title="Fidelity vs Success Probability\nfor different pulse durations")
	
	psa = Papprox.(ts,gᵒᵐ2,T₁ᵃ,γ2,nₜ2)
	infsa = 1 .- Fapprox.(ts,gᵒᵐ2,T₁ᵃ,γ2,nₜ2)
	m = infsa .< 1
	psa = psa[m]
	infsa = infsa[m]
	plot!(infsa,psa,ls=:dash,color=:black,label="simple approximation",
		legend=:bottomright,colorbar_tickfontsize=4)
end

# ╔═╡ a37026c0-a68a-4490-b408-1f54164d62c0
md"""# Final Summary

These results are used in the aforementioned publication. Consult the publication for discussion on how this informs the design of entanglement generation hardware and protocols.

There are interesting connections between temperature and mechanical quality factor discussed in the aforementioned publication and the rest of its supplementary materials."""

# ╔═╡ a3d54ce7-7f53-4cd4-a966-759d2fd9485e
begin
	#=
	# Old attempt at some additional plots
	import DataFrames, AlgebraOfGraphics, CairoMakie
	const df = DataFrames
	const ag = AlgebraOfGraphics
	const cm = CairoMakie
	
	records = []
	let
		npumps = range(1, 1000, length=10)
		Temps = range(0.1, 20, length=10) # K
		Qₘs = [1e3, 1e4, 1e5, 1e6]
		gᵒᵐ₀ = 0.268 # MHz
		ts = exp.(range(log(0.1),log(100),length=10)).*T₁ᵃ

		for npump in npumps, Temp in Temps, Qₘ in Qₘs, t in ts
			gᵒᵐ = sqrt.(npump) * gᵒᵐ₀ # MHz

			nₜ = n_per_T*Temp
			γ = π*F₀/Qₘ
			
			minT = Tminapprox(gᵒᵐ,T₁ᵃ,γ,nₜ)
	        if t < minT
				continue
			end
    		p = P(t,gᵒᵐ,T₁ᵃ,γ,nₜ)
			inf = 1 - F(t,gᵒᵐ,T₁ᵃ,γ,nₜ)
			
			push!(records, (;npump,Temp,Qₘ,gᵒᵐ,nₜ,γ,t,p,inf))
		end
	end
	dfrecords = df.DataFrame(records)
	
	summary = ag.data(dfrecords) * ag.mapping(:inf,:p) * ag.visual(cm.Lines)
	summary *= ag.mapping(color=:Temp, group=:Temp => ag.nonnumeric)
#	summary *= ag.mapping(row=:Qₘ => ag.nonnumeric, col)
#	summary *= ag.mapping(color=:Temp, group=:Qₘ => ag.nonnumeric)
	
	summary_plot = ag.draw(summary, axis=(xscale=log10,yscale=log10,cscale=log10))
	=#
end

# ╔═╡ 462de381-6e5e-4a92-bd1c-93a6fad5bfe0
let
	#=
	# Some rough estimates
	npump = 1:10.:1000.
	gᵒᵐ₀ = 0.268 # MHz
	gᵒᵐ = sqrt.(npump) * gᵒᵐ₀ # MHz
	
	Temp = 0.1:0.1:20 # K
	Qₘ2 = 10. ^ log₁₀Qₘ2
	nₜ2 = n_per_T*Temp2
	γ2 = π*F₀/Qₘ2

	md"""
	from Raniwala's design:
	$g_0$ = $(round(gᵒᵐ₀,digits=3)) MHz and T_a = $(round(T₁ᵃ*1000,digits=3)) ns
	(i.e. 1/Tₐ= $(round(1/T₁ᵃ/1000,digits=2)) GHz)
	
	derived:
	
	| n | g | T | Qₘ | nₜₕ |
	|:--|:--|:--|:---|:----|
	| $(round(n2)) | $(round(gᵒᵐ2,digits=1)) MHz | $(round(Temp2,digits=3)) K | $(round(Qₘ2,digits=-3)) | $(round(nₜ2,digits=2))
	"""
	=#
end

# ╔═╡ 74d20d02-3bb9-4cd5-ac08-143c68110e8d
md"""## Suggested Future Improvements

- What type of infidelity does this cause (what is the state obtained with probability 1-F)?
- This does not take into account the time necessary for cooling the mechanical resonator between attempts or the infidelity due to imperfect cooling.
- What happens when this is used to herald entangled pair (DLCZ)?
- What if the two nodes for the entanglement are not perfectly matched."""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuantumOptics = "6e0679c1-51ea-5a7c-ac74-d61b76210b0c"

[compat]
Plots = "~1.25.6"
PlutoUI = "~0.7.30"
QuantumOptics = "~1.0.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "91ca22c4b8437da89b030f08d71db55a379ce958"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.3"

[[Arpack_jll]]
deps = ["Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "e214a9b9bd1b4e1b4f15b22c0994862b66af7ff7"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.0+3"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "1ee88c4c76caa995a885dc2f22a5d548dfbbc0ba"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.2"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "bc1317f71de8dce26ea67fcdf7eccc0d0693b75b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.1"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CPUSummary]]
deps = ["Hwloc", "IfElse", "Static"]
git-tree-sha1 = "87b0c9c6ee0124d6c1f4ce8cb035dcaf9f90b803"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.6"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "926870acb6cbcf029396f2f2de030282b6bc1941"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.4"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "7b8f09d58294dc8aa13d91a8544b37c8a1dcbc06"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.4"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "6b6f04f93710c71550ec7e16b650c1b9a612d0b6"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.16.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DEDataArrays]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "31186e61936fbbccb41d809ad4338c9f7addf7ae"
uuid = "754358af-613d-5f8d-9788-280bf1605d4c"
version = "0.2.0"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DEDataArrays", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "15e43e11701b8c0b6250d7996b5768751f5a10c2"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.81.0"

[[DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "NLsolve", "OrdinaryDiffEq", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "e57ecaf9f7875714c164ccca3c802711589127cf"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.20.1"

[[DiffEqJump]]
deps = ["ArrayInterface", "Compat", "DataStructures", "DiffEqBase", "FunctionWrappers", "Graphs", "LinearAlgebra", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "628ddc7e2b44e214232e747b22f1a1d9a8f14467"
uuid = "c894b116-72e5-5b58-be3c-e6d8d4ac2b12"
version = "8.1.0"

[[DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "LinearAlgebra", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArrays", "Statistics"]
git-tree-sha1 = "d6839a44a268c69ef0ed927b22a6f43c8a4c2e73"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.9.0"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "9bc5dac3c8b6706b58ad5ce24cffd9861f07c94f"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.9.0"

[[Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "ab9e8f6b00c0584b103b7a8d4a221075a781845c"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.39"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExponentialUtilities]]
deps = ["ArrayInterface", "LinearAlgebra", "Printf", "Requires", "SparseArrays"]
git-tree-sha1 = "3e1289d9a6a54791c1ee60da0850f4fd71188da6"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.11.0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[FastBroadcast]]
deps = ["LinearAlgebra", "Polyester", "Static"]
git-tree-sha1 = "0f8ef5dcb040dbb9edd98b1763ac10882ee1ff03"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.12"

[[FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "6eae72e9943d8992d14359c32aed5f892bda1569"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.10.0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "2b72a5624e289ee18256111657663721d59c143e"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.24"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "4a740db447aae0fbeb3ee730de1afbb14ac798a1"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.63.1"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "aa22e1ee9e722f1da183eb33370df4c1aeb6c2cd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.63.1+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "d727758173afef0af878b29ac364a0eca299fc6b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.5.1"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HalfIntegers]]
git-tree-sha1 = "dc0ce9efc3d88c6cefc4e1f9c29b397be8734cfc"
uuid = "f0d1745a-41c9-11e9-1dd9-e5d34d218721"
version = "1.4.2"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8f0dc80088981ab55702b04bba38097a44a1a3a9"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.5"

[[Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "92d99146066c5c6888d5a3abc871e6a214388b91"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "2.0.0"

[[Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d8bccde6fc8300703673ef9e1383b11403ac1313"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.7.0+0"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "22df5b96feef82434b07327e2d3c770a9b21e023"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LRUCache]]
git-tree-sha1 = "d64a0aff6691612ab9fb0117b0995270871c5dfc"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.3.0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "41158dee1d434944570b02547d404e075da15690"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.7.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "83b56449c39342a47f3fcdb3bc782bd6d66e1d97"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.4"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LinearMaps]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "dbb14c604fc47aa4f2e19d0ebb7b6416f3cfa5f5"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.5.1"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SLEEFPirates", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "67c0dfeae307972b50009ce220aae5684ea852d1"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.101"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[ManualMemory]]
git-tree-sha1 = "9cb207b18148b2199db259adfa923b45593fe08e"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[NaNMath]]
git-tree-sha1 = "f755f36b19a5116bb580de457cda0c140153f283"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.6"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[NonlinearSolve]]
deps = ["ArrayInterface", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "200321809e94ba9eb70e7d7c3de8a7a6679a18b3"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.13"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "916077e0f0f8966eb0dc98a5c39921fdb8f49eb4"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.6.0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "Polyester", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "aab86872f5dccb1a6251ae7935b286e784dadba0"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "5.71.0"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ee26b350276c51697c9c2d88a072b339f9f03d73"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.5"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "92f91ba9e5941fc781fecf5494ac1da87bdac775"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.0"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "68604313ed59f0408313228ba09e79252e4b2da8"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.2"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "db7393a80d0e5bef70f2b518990835541917a544"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.6"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "5c0eb9099596090bb3215260ceca687b888a1575"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.30"

[[PoissonRandom]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "44d018211a56626288b5d3f8c6497d28c26dc850"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.0"

[[Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "a804f11a4da30c0561c824e09fc627531270d9a1"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.5.5"

[[PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "a3ff99bf561183ee20386aec98ab8f4a12dc724a"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.2"

[[PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "LabelledArrays"]
git-tree-sha1 = "e4cb8d4a2edf9b3804c1fb2c2de57d634ff3f36e"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.2.3"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[Primes]]
git-tree-sha1 = "984a3ee07d47d401e0b823b7d30546792439070a"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.1"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[QuantumOptics]]
deps = ["Arpack", "DiffEqCallbacks", "FFTW", "IterativeSolvers", "LinearAlgebra", "LinearMaps", "OrdinaryDiffEq", "QuantumOpticsBase", "Random", "RecursiveArrayTools", "Reexport", "SparseArrays", "StochasticDiffEq", "WignerSymbols"]
git-tree-sha1 = "eea423799266abf5c22a41313f26c11b9bb54488"
uuid = "6e0679c1-51ea-5a7c-ac74-d61b76210b0c"
version = "1.0.1"

[[QuantumOpticsBase]]
deps = ["Adapt", "FFTW", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "5506df67d351551b521c3099e19e06b6ea64e038"
uuid = "4f57444f-1401-5e15-980d-4471b28d5678"
version = "0.3.1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Random123]]
deps = ["Libdl", "Random", "RandomNumbers"]
git-tree-sha1 = "0e8b146557ad1c6deb1367655e052276690e71a3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.4.2"

[[RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[RationalRoots]]
git-tree-sha1 = "52315cf3098691c1416a356925027af5ab5bf548"
uuid = "308eb6b3-cc68-5ff3-9e97-c3c4da4fa681"
version = "0.2.0"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "37c1631cb3cc36a535105e6d5557864c82cd8c2b"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.0"

[[RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "6b96eb51a22af7e927d9618eaaf135a3520f8e2f"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.24.0"

[[RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "a6564a98066f512ff2efd438c8f1ce4262d69b87"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.7"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "62c2da6eb66de8bb88081d20528647140d4daa0e"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.0"

[[SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "1410aad1c6b35862573c01b96cd1f6dbe3979994"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.28"

[[SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "40c1c606543c0130cd3673f0dd9e11f2b5d76cd0"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.26.0"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "0afd9e6c623e379f593da01f20590bacc26d1d14"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.1"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "75c89362201983c500dd34923b015dbecdae7a90"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.20.0"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e08890d19787ec25029113e88c34ec20cac1c91e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.0.0"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2ae4fe21e97cd13efd857462c1869b73c9f61be3"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.2"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "d88665adc9bcf45903013af0982e2fd05ae3d0a6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "51383f2d367eb3b444c961d485c565e4c0cf4ba0"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.14"

[[StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "bedb3e17cc1d94ce0e6e66d3afa47157978ba404"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.14"

[[StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqJump", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "dd2043c1d182e11abf1b33e699dc3856e4663a54"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.42.0"

[[StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "Requires", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "12cf3253ebd8e2a3214ae171fbfe51e7e8d8ad28"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.2.9"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "d21f2c564b21a202f4677c0fba5b5ee431058544"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.4"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "884539ba8c4584a3a8173cb4ee7b61049955b79c"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.4.7"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "c3ab8b77b82fd92e2b6eea8a275a794d5a6e4011"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.9"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "Hwloc", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "6e261bff5c9f2537776165dea3067df9de4440cf"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.23"

[[VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "66d72dc6fcc86352f01676e8f0f698562e60510f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.23.0+0"

[[WignerSymbols]]
deps = ["HalfIntegers", "LRUCache", "Primes", "RationalRoots"]
git-tree-sha1 = "960e5f708871c1d9a28a7f1dbcaf4e0ee34ee960"
uuid = "9f57e263-0b3d-5e2e-b1be-24f2bb48858b"
version = "2.0.0"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─b576bacc-c769-4b0d-b033-ad1ad0cd4181
# ╟─74b628b2-82d5-471b-9f40-18b3b547668b
# ╟─89a7d7af-34a5-434d-b409-b40a070b5a3c
# ╟─a7d208cc-5c47-4fb3-8920-472e224bdc53
# ╟─e6b3962c-1a45-47ac-ad12-0e429c2bcccc
# ╟─aba7b5fd-9905-4915-a691-34cee7fbbb7b
# ╠═2fdaf940-174b-411f-905f-7d3dd101dc2f
# ╟─8f5f0424-1e03-473e-8cd4-32934fadc581
# ╠═5ea35653-e67d-4f33-9bbd-2f59d1c731a7
# ╟─05bcbadf-b764-4f0a-8115-8938860794ea
# ╟─a03f56e0-bf39-4645-a57f-0955903cf7a2
# ╟─69ad33a1-7497-4439-bdfa-8a8e30bd5068
# ╟─13430563-5293-4484-8b37-fcab0e442bd4
# ╟─87b86d27-520c-4947-91cb-ce321f245d26
# ╟─9c0a28b6-fb44-46cc-85d1-4c1daa4e098a
# ╠═72621af3-af5b-412d-86cb-a6376f983a5b
# ╠═5b6c3ce3-adfe-414d-b74f-65f4084f49ac
# ╟─35667cfb-b837-4197-b7d5-f59334109969
# ╟─4ff08761-cc25-4e04-9849-492332c39fd8
# ╠═78f9535c-9952-4b80-b681-73f9496371b1
# ╟─a37026c0-a68a-4490-b408-1f54164d62c0
# ╟─a3d54ce7-7f53-4cd4-a966-759d2fd9485e
# ╟─462de381-6e5e-4a92-bd1c-93a6fad5bfe0
# ╟─74d20d02-3bb9-4cd5-ac08-143c68110e8d
# ╟─637715ed-ceb6-4992-b665-026088559278
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
