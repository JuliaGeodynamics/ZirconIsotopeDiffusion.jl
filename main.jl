using CairoMakie
using NLsolve
using Interpolations
using Trapz
using CSV
using DataFrames
# using Parameters; export @kwdef, @unpack
using MuladdMacro
using BenchmarkTools, Test


abstract type NumericalParameters end
abstract type BoundaryParameters end
abstract type TraceElements end

@kwdef struct parameters{T,S,V} <:NumericalParameters
    Tsat::T           = 1032.5  # Starting at saturation
    Tend::T           = 695e0   # final temperature, C
    tfin::T           = 1500e0  # final ttime
    Cbulk::T          = 100e0
    RhoZrM::T         = 4.7/2.3 # Ratio of zircon to melt density
    Kmin::T           = 0.1     # Parameters for zirconium partition coefficient in major phase
    Kmax::T           = 0.1
    Crit::T           = 30e0
    delta::T          = 0.2
    Ktrace::T         = 0.1     # trace partition coefficient in major phase.
    Trace::S          = "Hf"
    XH20::T           = 2e0     # initial water content in melt, required for diffusion coefficient simulations.
    L::T              = 0.1     # 20e-4*(CZirc/CbulkZr)^(1./3.); radius of melt cel
    DGfZr::T          = 0.5     # diffusion coefficient ratio
    mass::V           = [89.9047026,90.9056439,91.9050386,93.9063148,95.908275]
end

@kwdef mutable struct BC_parameters{FP, V, S} <: BoundaryParameters
    D::V               = Float64[]
    csat::FP           = 0e0
    alpha::V           = Float64[]
    beta::V            = Float64[]
    T::FP              = 0e0
    Cz::FP             = 0e0
    Trace::S           = "Hf"
end

@kwdef struct Hf_TraceElement{_T} <: TraceElements
    a::_T = 11.29e3
    b::_T = 2.275
end

@kwdef struct Ti_TraceElement{_T} <: TraceElements
    a::_T = -11.05e3
    b::_T = 6.06
end

@kwdef struct Y_TraceElement{_T} <: TraceElements
    a::_T = 19.47
    b::_T = 13.04
end

@kwdef struct U_TraceElement{_T} <: TraceElements
    a::_T = 15.32
    b::_T = 9.17
end

@kwdef struct Th_TraceElement{_T} <: TraceElements
    a::_T = 13.02e3
    b::_T = 8.54
end

@kwdef struct Sm_TraceElement{_T} <: TraceElements
    a::_T = 13.338
    b::_T = -0.622
end

@kwdef struct Dy_TraceElement{_T} <: TraceElements
    a::_T = 2460.0
    b::_T = -0.867
end

@kwdef struct Yb_TraceElement{_T} <: TraceElements
    a::_T = 33460.0
    b::_T = -1.040
end

@kwdef struct P_TraceElement{_T} <: TraceElements
    a::_T = 7.646
    b::_T = 5.047
end


function KD_trace(T::_T, par::Hf_TraceElement{_T})  where {_T <: Real}
    return exp(par.a / T  - par.b) # true Kd_Hf from this model 2022
end #  5.804 ns (0 allocations: 0 bytes)

function KD_trace(T::_T, par::Ti_TraceElement{_T})  where {_T <: Real}
    return exp(par.a / T  + par.b) # true Kd_Hf from this model 2022
end

function KD_trace(T::_T, par::Y_TraceElement{_T})  where {_T <: Real}
    X = 1000.0 / T
    return exp(par.a * X - par.b)
end

function KD_trace(T::_T, par::U_TraceElement{_T})  where {_T <: Real}
    X = 1000.0 / T
    return exp(par.a * X - par.b)
end

function KD_trace(T::_T, par::Th_TraceElement{_T})  where {_T <: Real}
    return exp(par.a / T - par.b)
end

function KD_trace(T::_T, par::Sm_TraceElement{_T})  where {_T <: Real}
    Csat = ZrSaturation(T)
    KD = (par.a * Csat^(par.b))
    return KD
end

function KD_trace(T::_T, par::Dy_TraceElement{_T})  where {_T <: Real}
    Csat = ZrSaturation(T)
    KD = (par.a * Csat^(par.b))
    return KD
end #  29.332 ns (0 allocations: 0 bytes)

function KD_trace(T::_T, par::Yb_TraceElement{_T})  where {_T <: Real}
    Csat = ZrSaturation(T)
    KD = (par.a * Csat^(par.b))
    return KD
end

function KD_trace(T::_T, par::P_TraceElement{_T})  where {_T <: Real}
    X = 1000.0 / T
    return exp(par.a * X - par.b)
end

# Hf_trace = Hf_TraceElement()
# @btime KD_trace($1000.0,$Hf_trace) # 5.731
# Ti_trace = Ti_TraceElement()
# @btime KD_trace($1000.0,$Ti_trace) #5.807
# Y_trace = Y_TraceElement()
# @btime KD_trace($1000.0,$Y_trace) #6.363
# U_trace = U_TraceElement()
# @btime KD_trace($1000.0,$U_trace) # 6.118
# Th_trace = Th_TraceElement()
# @btime KD_trace($1000.0,$Th_trace) #5.734
# Sm_trace = Sm_TraceElement()
# @btime KD_trace($1000.0,$Sm_trace) #29.337
# Dy_trace = Dy_TraceElement()
# @btime KD_trace($1000.0,$Dy_trace) #29.316
# Yb_trace = Yb_TraceElement()
# @btime KD_trace($1000.0,$Yb_trace) # 28.672
# P_trace = P_TraceElement()
# @btime KD_trace($1000.0,$P_trace)  # 7.072



#Helper functions
function ZrSaturation(T::_T)  where {_T<: Real}# defining Zr saturation conditions
    # Csat = 4.414e7 / exp(13352/T) / 2 # Watson 96, Eq 1, in ppm Zr for checking. (divide by 1),or mol Zr (divide by 2)
    # Mfactor = 0.0000048*(T)^2 - 0.0083626*(T) + 4.8484463 # empirical relations from magma Fig.
    # differentiation calc (file  M_factorsforOleg.xlsx
    Mfactor=1.62;
    Csat=490000/exp(10108/T-1.16*(Mfactor-1)-1.48); # Boehnkeetal2013ChemGeol351,324 assuming 490000ppm in Zircon

    # Mfactor = 1.3
    # Csat = 490000 / exp(10108/T + 1.16*(Mfactor - 1) - 1.48) # Boehnkeetal2013ChemGeol351,324 assuming 490,000 ppm in Zircon
    # Csat = 490000 / (exp(12900/T - 0.85*(Mfactor - 1) - 3.80)) # Watson and Harrison 1983

    # for Monazite (Does not work for some reason):
    # H2O = 1wt% use below expression (Table 4, Rapp Watson 86)
    # Csat = 0.0000190 * exp(0.0143872 * T)
    # Csat = 600000 / (exp(-0.0144 * T + 24.177))
    # H2O = 6wt% use below expression (Table 4, Rapp Watson 86)
    # Csat = 0.00012 * exp(0.01352 * T)
    # Csat = 600000 / (exp(-0.0135 * T + 22.296))

    # for Apatite:
    # SiO2 = 0.68
    # Csat = 430000 / exp((-4800 + 26400 * SiO2) / T + 3.10 - 12.4 * SiO2) # Harrison Watson 84

    return Csat
end #1.484 ns (0 allocations: 0 bytes)

@btime ZrSaturation($(750 + 273.15))

function DiffusionCoefficient!(Di::AbstractArray, T::_T, x::_T, DGfZr::_T,  mass::AbstractArray) where {_T<:Real} # defining Zr diffusion coefficients in melts as f(T,X2O)
    theta = 1000.0 / T
    lnD   = -(11.4*x+3.13)/(0.84*x+1)-(21.4*x+47)/(1.06*x+1)*theta;  # best fit of Zr Diff coefficients (several workers) and WH83 dependence on XH2O
    # Dif   = exp(lnD)#*1e4*365*24*3600;  # in cm2/y
    cm2yr = 1e4*365*24*3600
    Dif   = cm2yr * exp(lnD)  # in cm2/y
    # bet = 0.05  # +0.059 Watkins et al., 2017; Méheut et al., 2021
    # Di = zeros(6)
    for i in 1:5
        Di[i]=Dif*(mass[1]/mass[i])^0.05
    end
    # lnHf = -3.52 - 231.09 / 8.31 / theta
    # lnD_Hf = (-8.620340372 * T - 42705.17449 - .318918919 * x * T + 4049.500765 * x) / T
    Di[6] = Dif[1] * DGfZr  # exp(lnD_Hf) * 1e4 * 365 * 24 * 3600  # in cm2/y
    return nothing
end   #112.346 ns (2 allocations: 32 bytes)

Di = zeros(6)
parameter_test = parameters()
@btime DiffusionCoefficient!($Di, $1023.15, $2.0, $0.5, $parameter_test.mass)
DiffusionCoefficient!(Di, 1023.15, 2.0, 0.5, parameter_test.mass)
@test sum(Di) ≈ 6.6006e-5 atol = 1e-5

function bc(X::_T, Eq::_T, BC_parameters::BoundaryParameters, TraceElement::TraceElements) where {_T}
    ct = BC_parameters.alpha .* X +BC_parameters.beta
    grad = -BC_parameters.D .* (ct - X)
    # Eq = zeros(6)
    Eq[1] = sum(X[i] for i in 1:5) - BC_parameters.csat
    @. Eq[2:5] = grad[2:5] * X[1] - X[2:5] * grad[1]
    # @. Eq[2:5] = grad[2:5] * X[1] - X[2:5]' * grad[1] #matlab version
    KD_Hf = KD_trace(BC_parameters.T, TraceElement)
    # KD_Hf = kdHf(BC_parameters.T, BC_parameters)
    CHfs = X[6] * KD_Hf
    Cz = BC_parameters.Cz * X[1] / BC_parameters.csat
    Eq[6] = Cz * grad[6] - grad[1] * (CHfs - X[6])
    return Eq
end #2.548 ns (0 allocations: 0 bytes)

Eq = zeros(6)
BC_parameter_test = BC_parameters()
TraceElement_test = Hf_TraceElement()
@btime f = (X) -> bc(X, Eq, $BC_parameter_test, $TraceElement_test)

function mf_magma(Tk::_T) where {_T <: Real}
    T = Tk - 273.15
    t2 = T * T
    t7 = exp(@muladd 0.961026371384066e3 - 0.3590508961e1 * T + (0.4479483398e-2- 0.1866187556e-5 * T) * t2)
    CF = inv(1.0 + t7)
    return CF
end #18.418 ns (0 allocations: 0 bytes)

@btime mf_magma($1000)
@test mf_magma(1000.0) ≈ 0.23081 atol = 1e-4

function mf_rock(T::_T) where {_T <: Real}
    t2 = T * T
    t7 = exp(@muladd 0.961026371384066e3 - 0.3590508961e1 * T + (0.4479483398e-2- 0.1866187556e-5 * T) * t2)
    CF = inv(1.0 + t7)
    return CF
end #20.516 ns (0 allocations: 0 bytes)

@btime mf_rock($1000)
@test mf_rock(1000.0) ≈ 0.99999 atol = 1e-4

function progonka(C0, dt, it, Di, Xs, Temp, MeltFrac, Dplag, Dscale, Dflux, parameter::NumericalParameters, BC_parameters::BoundaryParameters, TraceElement::TraceElements; A=A, B=B, C=C, D=D, F=F, alpha=alpha, beta= beta, Eq=Eq, n=500, R = range(0, stop=1, length=n))
    CZirc = BC_parameters.Cz;
    Temp = Temp[it];
    MeltFrac_new = MeltFrac[it+1];
    MeltFrac_old = MeltFrac[it-1];
    MeltFrac = MeltFrac[it];
    S = (Xs^3 + MeltFrac*(1.0 - Xs^3))^(1.0/3.0); # rad of the melt shell

    DiffusionCoefficient!(Di, Temp, parameter.XH20, parameter.DGfZr, parameter.mass)
    Dif = Di
    Dif ./=Dscale; #see below Diff Coeff dependednt on water and T in cm2/s
    Csat = ZrSaturation(Temp);

    Czl = CZirc*C0[1,1]/Csat;
    Czh = CZirc*C0[1,4]/Csat;

    cte = inv(R[2]-R[1])
    ctf = inv(S-Xs)
    V = 0.0
    for j in 1:5
        @inbounds Dflux[j] = a = Dif[j]*(C0[2,j] - C0[1,j]) * cte * ctf;
        V += a
    end
    V = -V/(CZirc*parameter.RhoZrM-Csat);
    # @. Dflux[1:5]=Dif[1:5]*(C0[2,1:5] - C0[1,1:5])* cte *ctf;
    # V=-sum(Dflux)/(CZirc*parameter.RhoZrM-Csat);

    diffF= if it>1
        (MeltFrac_new-MeltFrac_old)*inv(dt)*0.5;
    else
        (MeltFrac_new-MeltFrac)*inv(dt);
    end

    @muladd W=(1/3)*(diffF*(1-Xs^3)-3*Xs^2*V*(MeltFrac-1))/((-MeltFrac+1)*Xs^3+MeltFrac)^(2/3);
    dC=sum(C0[n,i] for i in 1:5)-Csat;
    t4 = tanh(parameter.delta * (dC - parameter.Crit));
    t7 = tanh(parameter.delta * parameter.Crit);
    cte =  @muladd 0.1e1 / (0.1e1 + t7) * (t4 * (parameter.Kmax - parameter.Kmin) + parameter.Kmax * t7 + parameter.Kmin);
    for j in 1:5
        Dplag[j] = cte
    end
    Dplag[6]=parameter.Ktrace;

    for j in axes(D,2)
        Dif_j = Dif[j]
        D[n,j] = @muladd -Dif_j-W*(R[end]-R[end-1])*(S-Xs)*(1-Dplag[j]);
        A[n,j] = Dif_j
        F[n,j] = 0;
    end
    # Coefficients for Thomas method
    s = Xs
    for j in axes(F,2)
        @inbounds for i in 2:n-1
            psi1 = R[i-1]
            psi2 = R[i]
            psi3 = R[i+1]
            t1 = Dif[j] * dt
            @muladd t5 = (psi1 * S - psi1 * s + s) ^ 2
            t6 = psi2 - psi1
            t8 = t5 / t6
            t12 = S * psi2
            @muladd t14 = ((-psi2 + 1) * s + t12) ^ 2
            t15 = S - s
            t20 = (-W + V) * psi2 - V
            @muladd A[i,j] = -t14 * t15 * dt * psi2 * t20 - t1 * t8
            @muladd t25 = (-psi2 * s + s + t12) ^ 2
            t28 = t25 / (psi3 - psi2)
            B[i,j] = -t1 * t28
            t32 = -t15
            t33 = t32 ^ 2
            t34 = -t6
            @muladd t38 = (t32 * psi2 - s) ^ 2
            @muladd D[i,j] = -t1 * (-t28 - t8) - t33 * t34 * t38 - t20 * psi2 * dt * t38 * t32
            t44 = t15 ^ 2
            @muladd t48 = (t15 * psi2 + s) ^ 2
            @muladd F[i,j] = -t34 * t44 * t48 * C0[i,j]
        end
    end

    # Forward Thomas path
    @inbounds for j in axes(alpha, 2)
        c = inv(D[n, j])
        alpha[n, j] = -A[n, j] * c
        beta[n, j]  =  F[n, j] * c
    end

    @inbounds for j in axes(alpha, 2), i in n-1:-1:2
        Bij = B[i,j]
        Dij = D[i,j]
        alphaij = alpha[i+1,j]
        alpha[i, j] = @muladd -A[i, j] / (Bij * alphaij + Dij)
        beta[i, j]  = @muladd (F[i,j] - Bij * beta[i+1,j]) / (Bij * alphaij + Dij)
    end

    # Boundary conditions
    BC_parameters.D = Dif[:]
    BC_parameters.csat = Csat
    BC_parameters.alpha = alpha[2,:]
    BC_parameters.beta = beta[2,:]
    BC_parameters.T = Temp
    # BC_parameters.Cz = CZirc
    # @show CZirc
    # BC_parameters.Trace = parameter.Trace

    f = (X) -> bc(X, Eq, BC_parameters, TraceElement) # function of dummy variable y
    result = NLsolve.nlsolve(f, (@views C0[1,:]), method = :trust_region) #NLsolve doesnt provide the Levenberg-Marquart method, but trust_region comes close to it
    out = result.zero  # solution vector

    # if result.f_converged <= 0  #convergence
    #     println(result.residual_norm) #residual norm
    # end
    @views C[1,:] = out[:]

    # Backward Thomas path
    for j in axes(C,2), i in 1:n-1
        @inbounds C[i+1,j] = @muladd C[i,j] * alpha[i+1,j] + beta[i+1,j]
    end

    return C, Czl, Czh, Csat, Dif, S, Dplag, V, W

end

function TemperatureHistory_m_erupt(tr, Tr, nt, par::NumericalParameters)
    if isempty(tr)
        ti = range(0, stop=par.tfin, length=nt) |> collect
        Ti = range(par.Tsat, stop=par.Tend+273.15, length=nt) |> collect
        # CrFrac1 = mf_rock.(Ti .- 273.15)
        CrFrac1 = [mf_rock(Ti - 273.15) for Ti in Ti]
    else
        istart = findfirst(x -> x > 0, Tr)
        Tr[istart-1] = 950 + 273.15
        tr[istart-1] = tr[istart] - 5
        dT = 0.05
        if minimum(mf_rock(Ti) for Ti in Tr) < 0.01
        # if minimum(mf_rock.(Tr)) < 0.01
            println("no melting")
            return Float64[], Float64[], Float64[]
        end
        if minimum(Tr) - Tsat > 0
            println("high temperature")
            return Float64[], Float64[], Float64[]
        end
        ttime = tr
        Temp = Tr
        try
            it = findfirst(x -> x < Tsat, Temp)
            ttime[it-1] = ttime[it] - (Temp[it] - Tsat) / (Temp[it] - Temp[it-1]) * 5
            Temp[it-1] = Tsat
            time1 = ttime[it-1:end]
            Temp1 = Temp[it-1:end]
            nt = length(time1)
            s = zeros(Temp1)
            for i in 2:nt
                s[i] = s[i-1] + abs(Temp1[i] - Temp1[i-1])
            end
            ni = floor(s[nt] / dT)
            si = range(s[1], stop=s[nt], length=ni)
            ti = interp1(s, time1, si)
            Ti = interp1(time1, Temp1, ti)
        catch ME
            # println("wrong Thist for sample: ", sampnum, ", ", ME.message)
            return Float64[], Float64[], Float64[]
        end
        CrFrac1 = mf_rock.(Ti .- 273.15)
    end
    return ti, Ti, CrFrac1
end

tr_test = Float64[]
Tr_test = Float64[]
@btime TemperatureHistory_m_erupt($tr_test, $Tr_test,$500, $parameter_test)

@views function ZirconIsotopeDiffusion(; n = 500, nt = 500, tyear = (3600*24*365))

    Runname = "Test"
    !isdir("Results") && mkpath("Results")

    # parameters for simulations
    CbulkZr = 100.0
    tyear = 3600*24*365
    iplot = 1 # plot results
    n = 500 # number of points in radial mesh. Can be changed by user depends on desired accuracy
    nt = 500
    CZirc = 490000.0 # zirconium concentration in zircon, ppm
    XH2O = 2.0 # initial water content in melt, required for diffusion coefficient simulations.
    Tsolidus = 400 + 273 # arbitrary solidus temperature for phase diagram used
    Csupsat = 3 # ppm supersaturation to cause nucleation of a new crystal upon cooling
    UCR = 1 # Critical concentration for microZircon nucleation on major minerals
    ZircNuc = 1e-4 # Zircon stable nuclei in cm
    L = 0.1 # 20e-4*(CZirc/CbulkZr)^(1./3.); radius of melt cell
    DGfZr = 0.5 # ratio of diffusion coefficients of Hf to Zr; change for other element of interest
    # mass=[89.9047026,90.9056439,91.9050386,93.9063148,95.908275]; #masses defined by O.Melnik's code

    # Allocations (formerly matrixes function)
    C0 = zeros(n, 6)
    C = zeros(n, 6)
    A = zeros(n, 6)
    B = zeros(n, 6)
    D = zeros(n, 6)
    F = zeros(n, 6)
    alpha = zeros(n, 6)
    beta = zeros(n, 6)
    x = range(0, 1, n)
    VV = zeros(nt, 1)  # arrays for future storage of data and plotting
    XXs = zeros(nt, 1)
    RRd = zeros(nt, 1)
    tt = zeros(nt, 1)
    UU = zeros(nt, 1)  # array for undersaturation from first to last distance length point
    Tsave = zeros(nt, 1)
    ZrPls = zeros(nt, 1)
    Xp_sav = zeros(nt, 1)
    CC = zeros(nt, n, 6);
    Dplag = zeros(1,6)
    Zcomp = zeros(1, nt)
    ZrHF = zeros(1, nt)
    CZircon = zeros(1, 5)
    Cplag = zeros(1, 5)
    CintS = zeros(n-1, 5)
    Cint = zeros(1,5)
    Dflux = zeros(1, 5)
    Zcompl = zeros(1, nt-1)
    Zcomph = zeros(1, nt-1)
    Melndelta = zeros(1, nt-1)
    Di = zeros(6)
    Eq = zeros(6)
    R = range(0, stop=1, length=n)
    rr = range(0, stop=1, length=n)


    # Solve for Tsat
    function equation!(F, T)
        F[1] = ZrSaturation(T[1])*mf_rock.(T[1]-273.15) - CbulkZr
    end

    result = nlsolve(equation!, [1000.0])
    Tsat = result.zero[1]

    # parameters for the simulation
    parameter = parameters(;
        Tsat = Tsat, # Starting at saturation
        Tend = 695e0,  # final temperature, C
        tfin = 1500e0, # final time
        Cbulk = CbulkZr,
        RhoZrM = 4.7/2.3, # Ratio of zircon to melt density
        Kmin = 0.1,  # Parameters for zirconiun partition coefficient in major phase
        Kmax = 0.1,
        Crit = 30e0,
        delta = 0.2,
        Ktrace = 0.1, # trace partition coefficient in major phase.
        Trace = "Hf",
        XH20 = XH2O,
        L = L,
        DGfZr = DGfZr,  # diffusion coefficient ratio
    )
    tr = Float64[]
    Tr = Float64[]
    ttime, Temp, MeltFrac = TemperatureHistory_m_erupt(tr, Tr, nt, parameter)

    # BC parameters
    BC_parameter = BC_parameters(;
    Cz = CZirc,
    Trace = parameter.Trace,
    )
    TraceElement = Hf_TraceElement()
    # Scaling
    tfin = parameter.tfin[end] # total time in years of the process
    # SCALING-----------------
    DiffusionCoefficient!(Di, (750 + 273.15), XH2O, DGfZr, parameter.mass)
    Ds = Di
    Dscale = Ds[1]
    tscale = L^2 / Dscale # dimensionless scale for the time
    ttime = ttime ./ tscale
    # nt = L(ttime)  # this is obsolete as the nt does not change with scaling
    # END:SCALING-----------------

    # Initial Conditions
    t = ttime[1] / tscale
    ZirconRadius = 2e-4
    Xs = ZirconRadius / L
    ZircNuc = ZircNuc / L
    S0 = S = (Xs^3 + MeltFrac[1] * (1 - Xs^3))^(1/3)
    dt = ttime[2] - ttime[1]
    W = 0
    V = 0


    @views C0[:, 1] .= ZrSaturation(Temp[1]) * 0.5145
    @views C0[:, 2] .= ZrSaturation(Temp[1]) * 0.1122
    @views C0[:, 3] .= ZrSaturation(Temp[1]) * 0.1715
    @views C0[:, 4] .= ZrSaturation(Temp[1]) * 0.1738
    @views C0[:, 5] .= ZrSaturation(Temp[1]) * 0.0280
    @views C0[:, 6] .= CZirc / KD_trace(Temp[1], TraceElement) / 70
    # C0[1:n,6] = 50  # PHOSHPORUS< CHANGEHF melt from Bachmann etal JPet 2002.
    Dplag[1:5] .= 0.1
    Dplag[6] = 0.1
    sleep(1e-5)

    CC[1, 1:n, 1:5] = C0[1:n, 1:5]
    Tsave[1] = Temp[1] - 273.15
    XXs[1] = Xs * 1e4 * L
    RRd[1] = S * 1e4 * L
    ZrPls[1] = XXs[1] # zircon radius in um
    UU[1] = C0[1, 1]
    tt[1] = ttime[1] * tscale
    Zcomp[1,1] = C0[1, 4] / C0[1, 1]
    ZrHF[1,1] = CZirc / KD_trace(Temp[1], TraceElement) / C0[1, 6]
    # Zcomp[1] = C0[1, 4] / C0[1, 1] #matlab version
    # ZrHF[1] = CZirc / KD_trace(Temp[1], par) / C0[1, 6] #matlab version
    Melndelta[1,1] = Zcomp[1,1]


    @views @. CZircon[1:5] = 4 * π * CZirc * C0[1, 1:5] / ZrSaturation(Temp[1]) * ZirconRadius^3 / 3
    @views Cplag[1:5] .= 0
    @views @. CintS[1, 1:5] = CZircon[1:5] + 4 * π * C0[1, 1:5] * (S^3 - ZirconRadius^3) / 3

    # Main loop
    @views for i = 2:nt-1
        if MeltFrac[i] > 0.01
            C, Czl, Czh, Csat, Dif, S, Dplag, V, W = progonka(C0,dt,i,Di,Xs,Temp,MeltFrac, Dplag, Dscale,Dflux,parameter, BC_parameter, TraceElement; A=A, B=B, C=C, D=D,F=F, alpha=alpha, beta=beta,Eq=Eq, n=n, R=R)
            dt = ttime[i] - ttime[i-1]
            C0 = C
        else
            V = 0
            W = 0
        end
        rr = R * (S - Xs) .+ Xs
        Csat = ZrSaturation(Temp[i])
        @views CZircon[1:5] = CZircon[1:5] - CZirc * C[1, 1:5] / Csat * 4 * π * Xs^2 * V * dt
        @views Cplag[1:5] = Cplag[1:5] - C[end, 1:5] .* Dplag[1:5] * 4 * π * S^2 * W * dt
        @views Cint[1:5] .= 0
        for ik = 2:n
            for i in 1:5
                Cint[i] = Cint[i] + (C[ik-1, i] * rr[ik-1]^2 + C[ik, i] * rr[ik]^2) / 2 * (rr[ik] - rr[ik-1])
            end
        end
        Cint = 4 * π * Cint + parameter.RhoZrM * CZircon + Cplag
        @views CintS[i, 1:5] = Cint[1:5]

        if iplot == 1 && i % floor(nt / 10) == 0
            fig = Figure(size = (800, 800), backgroundcolor = :white)

            # Subplots
            ax1 = Axis(fig[1, 1], xlabel = "Distance, um", ylabel = L"\delta^{94/90}Zr")
            ax2 = Axis(fig[2, 1], xlabel = "Distance )", ylabel = "Zr/Hf")

            rr = R * (S - Xs) .+ Xs

            # Plot data (replace `data` with your actual data)
            lines!(ax1, rr * L * 1e4, (C[:, 4] .* 0.5145 ./ C[:, 1] ./ 0.1738 .- 1) .* 1000, linewidth = 1.5)
            lines!(ax2, rr * L * 1e4, (sum(C[:, 1:5], dims = 2) ./ C[:, 6])[:], linewidth = 1.5)
            display(fig)
        end


        t += dt
        rr = R * (S - Xs) .+ Xs
        Cl = trapz(rr[:,1], rr[:,1].^2 .* C[:, 1])
        Ch = trapz(rr, rr.^2 .* C[:, 4])
        Xs = max(ZircNuc, Xs - V * dt)
        S0 = S
        XXs[i] = Xs * 1e4 * L  # zircon radius in um
        RRd[i] = S * 1e4 * L  # melt cell radius in um
        VV[i] = -V * L * 1e4 / tscale  # array of dissolution rate
        tt[i] = ttime[i] * tscale
        UU[i] = C[1] - ZrSaturation(Temp[i])
        Tsave[i] = Temp[i] - 273
        ZrPls[i] = minimum(XXs[1:i, 1])
        Zcompl[i] = Czl / CZirc
        Zcomph[i] = Czh / CZirc
        Zcomp[1,i] = C[1, 4] / C[1, 1]
        Melndelta[1,i] = Ch / Cl
        ZrHF[i] = CZirc / KD_trace(Temp[i], TraceElement) / C0[1, 6]
        @views CC[i, 1:n, 1:6] .= C0[1:n, 1:6]
    end
    # Plot results (if iplot is set)



    if iplot == 1
        fig = Figure(size = (800, 800), backgroundcolor = :white)

        # Subplots
        ax1 = Axis(fig[1, 1], xlabel = "Time (years)", ylabel = "Zr radius")
        ax2 = Axis(fig[1, 2], xlabel = "Time (years)", ylabel = "Temperature T, ^oC")
        ax3 = Axis(fig[2, 1], xlabel = "Distance", ylabel = L"Growth Rate, cm.s^{-1}")
        ax4 = Axis(fig[2, 2], xlabel = "Distance, um", ylabel = L"\delta^{94/90}Zr")
        ax5 = Axis(fig[3, 1], xlabel = "Distance ", ylabel = "Zr/Hf")

        # Plot data (replace `data` with your actual data)
        lines!(ax1, tt[1:end-1]/1e3, XXs[1:end-1], color = :blue)
        lines!(ax2, tt[1:end-1]/1e3, Tsave[1:end-1], color = :blue)
        lines!(ax3, XXs[1:end-1], VV[1:end-1], color = :blue)
        DelZr = zeros(1,nt)
        DelZr[2:end-1] = (Zcomp[2:end-1] ./ Zcomp[2] .- 1) * 1000
        DelMlt = (Melndelta[2:end-1] ./ Melndelta[2] .- 1) * 1000
        lines!(ax4, XXs[1:end-1], DelZr[1:end-1], color = :blue)
        lines!(ax5, XXs[2:end-1], ZrHF[2:end-1], color = :blue)

        display(fig)
    end

    # # Save results

    # # Print the figure to a PDF file
    # CairoMakie.save("Results/Test.pdf", fig)

    # i = nt - 2

    # # Convert array to DataFrame
    # Rsave = DataFrame(time_ka = tt[1:i-1] / 1e3, Rad_um = XXs[1:i-1], Gr_rate_mm_a = VV[1:i-1], Temp_C = Tsave[1:i-1], DelZr = DelZr[1:i-1], DelMlt = DelMlt[1:i-1], ZrHf = ZrHF[1:i-1])

    # # Write DataFrame to CSV file
    # CSV.write("Results/$Runname.csv", Rsave)

    # # Append structure to CSV file
    # par = DataFrame(fname = "$Runname")

    # if !isfile("Results/summary.csv")
    #     wwar = true
    # else
    #     wwar = false
    # end

    # if wwar
    #     CSV.write("Results/summary.csv", par)
    # else
    #     existing = CSV.read("Results/summary.csv", DataFrame)
    #     append!(existing, par)
    #     CSV.write("Results/summary.csv", existing)
    # end

    # return C, XXs, DelZr
end

@code_warntype ZirconIsotopeDiffusion()
@btime ZirconIsotopeDiffusion()
ZirconIsotopeDiffusion()

# #Matlab
# # >> sum(C(1,:))
# # ans = 129.4123
# @test sum(C[1,:]) ≈ 129.4123 atol=1e-4

# # >> max(XXs)
# # ans = 19.4589
# @test maximum(XXs) ≈ 19.4589 atol=1e-4

# # >> min(DelZr)
# # ans = -2.0528
# # Julia = -2.5351 (tolerance??)
# @test minimum(DelZr[2:end-1]) ≈ -2.0528 atol=1e-3
