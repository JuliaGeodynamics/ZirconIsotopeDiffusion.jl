using NLsolve
using CairoMakie
using Interpolations
using Trapz
using CSV
using DataFrames
using Test
using BenchmarkTools

#Helper functions
function ZrSaturation(T) # defining Zr saturation conditions
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
end

function DiffusionCoefficient(T, x, DGfZr)  # defining Zr diffusion coefficients in melts as f(T,X2O)
    # global DGfZr
    theta = 1000 / T
    lnD=-(11.4*x+3.13)/(0.84*x+1)-(21.4*x+47)/(1.06*x+1)*theta;  # best fit of Zr Diff coefficients (several workers) and WH83 dependence on XH2O
    Dif=exp(lnD)*1e4*365*24*3600;  # in cm2/y
    mass=[89.9047026,90.9056439,91.9050386,93.9063148,95.908275];
    bet = 0.05  # +0.059
    Di = zeros(6)
    Di[1:5]=Dif*(mass[1]./mass).^bet;
    # lnHf = -3.52 - 231.09 / 8.31 / theta
    # lnD_Hf = (-8.620340372 * T - 42705.17449 - .318918919 * x * T + 4049.500765 * x) / T
    Di[6] = Dif[1] * DGfZr  # exp(lnD_Hf) * 1e4 * 365 * 24 * 3600  # in cm2/y
    return Di
end

@test sum(DiffusionCoefficient((750 + 273.15), 2, (0.5))) ≈ 6.60069e-5 atol=1e-4

function kdHf(T, par)
    X = 1000. / T
    KD_Hf = exp(11.29e3 / T  - 2.275) # true Kd_Hf from this model 2022
    KD_Ti = exp(-11.05e3 / T + 6.06) # Kd for Ti zircon/melt based on Ferry and Watson 2007
    KD_Y = exp(19.47 * X - 13.04)
    KD_U = exp(15.32 * X - 9.17) #U
    KD_Th = exp(13.02e3 / T -8.54) # Kd for Th
    Csat = ZrSaturation(T)
    KD_Sm = (13.338 * Csat^(-0.622))
    KD_Dy = (2460.0 * Csat^(-0.867))
    KD_Yb = (33460. * Csat^(-1.040))
    KD_P = exp(7.646 * X - 5.047)

    KD = get(
        Dict(
            "Hf" => KD_Hf,
            "Y" => KD_Y,
            "U" => KD_U,
            "P" => KD_P,
            "Sm" => KD_Sm,
            "Dy" => KD_Dy,
            "Yb" => KD_Yb,
            "Th" => KD_Th,
            "Ti" => KD_Ti
        ),
        par["Trace"],
        nothing
    )

    return KD
end

@test kdHf(1023.15, Dict("Trace" => "Hf")) ≈ 6371.245 atol=1e-3
@btime kdHf($1023.15, $Dict("Trace" => "Hf"))

function bc(X, par)
    ct = par["alpha"] .* X + par["beta"]
    grad = -par["D"] .* (ct - X)
    Eq = zeros(6)
    Eq[1] = sum(X[1:5]) - par["csat"]
    @. Eq[2:5] = grad[2:5] * X[1] - X[2:5] * grad[1]
    # @. Eq[2:5] = grad[2:5] * X[1] - X[2:5]' * grad[1] #matlab version
    KD_Hf = kdHf(par["T"], par)
    CHfs = X[6] * KD_Hf
    Cz = par["Cz"] * X[1] / par["csat"]
    Eq[6] = Cz * grad[6] - grad[1] * (CHfs - X[6])
    return Eq
end

function mf_magma(Tk)
    T = Tk - 273.15
    t2 = T .* T
    t7 = exp.(0.961026371384066e3 .- 0.3590508961e1 .* T .+ 0.4479483398e-2 .* t2 .- 0.1866187556e-5 .* t2 .* T)
    CF = 0.1e1 ./ (0.1e1 .+ t7)
    return CF
end

@btime mf_magma($1000)
@test mf_magma(1000) ≈ 0.23081 atol = 1e-4

function mf_rock(T)
    t2 = T .* T
    t7 = exp.(0.961026371384066e3 .- 0.3590508961e1 .* T .+ 0.4479483398e-2 .* t2 .- 0.1866187556e-5 .* t2 .* T)
    CF = 0.1e1 ./ (0.1e1 .+ t7)
    return CF
end

@btime mf_rock($1000)
@test mf_rock(1000) ≈ 0.99999 atol = 1e-4

function progonka(C0,dt,it,parameters)

    global  n, R,Dplag, ZrPl, MinCore, time, tscale, S0
    global A, B, D, F, alpha, beta, Xs, Temp, MeltFrac, XH2O, Tsolidus, V, W, Csupsat, Dscale, UCR, CZirc, S, ZircNuc, Czl, Czh, Dflux

    S=(Xs^3+MeltFrac[it]*(1-Xs^3))^(1/3); # rad of the melt shell
    Dif=DiffusionCoefficient(Temp[it],XH2O, DGfZr)/Dscale; #see below Diff Coeff dependednt on water and T in cm2/s
    Csat=ZrSaturation(Temp[it]);
    Czl=CZirc*C0[1,1]/Csat;
    Czh=CZirc*C0[1,4]/Csat;

    Dflux[1:5]=Dif[1:5].*(C0[2,1:5] - C0[1,1:5])/(R[2]-R[1])/(S-Xs);
    V=-sum(Dflux)/(CZirc*parameters["RhoZrM"]-Csat);

    if it>1
        diffF=(MeltFrac[it+1]-MeltFrac[it-1])/dt/2;
    else
        diffF=(MeltFrac[it+1]-MeltFrac[it])/dt;
    end

    W=(1/3)*(diffF*(1-Xs^3)-3*Xs^2*V*(MeltFrac[it]-1))/((-MeltFrac[it]+1)*Xs^3+MeltFrac[it])^(2/3);
    dC=sum(C0[n,1:5])-Csat;
    Ccr=parameters["Crit"];
    delta=parameters["delta"];
    Dpmax=parameters["Kmax"];
    Dpmin=parameters["Kmin"];
    t4 = tanh(delta * (dC - Ccr));
    t7 = tanh(delta * Ccr);
    @. Dplag[1:5] = 0.1e1 / (0.1e1 + t7) * (t4 * (Dpmax - Dpmin) + Dpmax * t7 + Dpmin);
    Dplag[6]=parameters["Ktrace"];

    @. D[n,:]=-Dif[:]-W*(R[n]-R[n-1])*(S-Xs)*(1-Dplag[:]);
    @. A[n,:]=Dif[:];
    @. F[n,:]=0;
    # Coefficients for Thomas method
    s = Xs
    for j in 1:6
        for i in 2:n-1
            psi1 = R[i-1]
            psi2 = R[i]
            psi3 = R[i+1]
            t1 = Dif[j] * dt
            t5 = (psi1 * S - psi1 * s + s) ^ 2
            t6 = psi2 - psi1
            t8 = t5 / t6
            t12 = S * psi2
            t14 = ((-psi2 + 1) * s + t12) ^ 2
            t15 = S - s
            t20 = (-W + V) * psi2 - V
            A[i,j] = -t14 * t15 * dt * psi2 * t20 - t1 * t8
            t25 = (-psi2 * s + s + t12) ^ 2
            t28 = t25 / (psi3 - psi2)
            B[i,j] = -t1 * t28
            t32 = -t15
            t33 = t32 ^ 2
            t34 = -t6
            t38 = (t32 * psi2 - s) ^ 2
            D[i,j] = -t1 * (-t28 - t8) - t33 * t34 * t38 - t20 * psi2 * dt * t38 * t32
            t44 = t15 ^ 2
            t48 = (t15 * psi2 + s) ^ 2
            F[i,j] = -t34 * t44 * t48 * C0[i,j]
        end
    end

    # Forward Thomas path
    alpha[n,:] = -A[n,:] ./ D[n,:]
    beta[n,:] = F[n,:] ./ D[n,:]
    for i in n-1:-1:2
        alpha[i,:] = -A[i,:] ./ (B[i,:] .* alpha[i+1,:] + D[i,:])
        beta[i,:] = (F[i,:] - B[i,:] .* beta[i+1,:]) ./ (B[i,:] .* alpha[i+1,:] + D[i,:])
    end

    # Boundary conditions
    parb = Dict()
    parb["D"] = Dif[:]
    parb["csat"] = Csat
    parb["alpha"] = alpha[2,:]
    parb["beta"] = beta[2,:]
    parb["T"] = Temp[it]
    parb["Cz"] = CZirc
    parb["Trace"] = parameters["Trace"]

    f = (X) -> bc(X, parb) # function of dummy variable y
    result = NLsolve.nlsolve(f, C0[1,:], method = :trust_region) #NLsolve doesnt provide the Levenberg-Marquart method, but trust_region comes close to it
    out = result.zero  # solution vector
    fval = result.residual_norm  # residual vector
    exflag = result.f_converged  # convergence flag (true if converged)

    if exflag <= 0
        println(fval)
    end
    C[1,:] = out[:]

    # Backward Thomas path
    for i in 1:n-1
        C[i+1,:] = C[i,:] .* alpha[i+1,:] + beta[i+1,:]
    end

    return C, Czl, Czh, Csat, Dif, S

end

function TemperatureHistory_m_erupt(tr, Tr, par)
    if isempty(tr)
        nt = par["nt"]
        ti = range(0, stop=par["tfin"], length=nt)
        Ti = range(par["Tsat"], stop=par["Tend"]+273.15, length=nt)
        CrFrac1 = mf_rock(Ti .- 273.15)
    else
        istart = findfirst(x -> x > 0, Tr)
        Tr[istart-1] = 950 + 273.15
        tr[istart-1] = tr[istart] - 5
        dT = 0.05
        if minimum(mf_rock(Tr)) < 0.01
            println("no melting")
            return [], [], []
        end
        if minimum(Tr) - Tsat > 0
            println("high temperature")
            return [], [], []
        end
        time = tr
        Temp = Tr
        try
            it = findfirst(x -> x < Tsat, Temp)
            time[it-1] = time[it] - (Temp[it] - Tsat) / (Temp[it] - Temp[it-1]) * 5
            Temp[it-1] = Tsat
            time1 = time[it-1:end]
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
            println("wrong Thist for sample: ", sampnum, ", ", ME.message)
            return [], [], []
        end
        CrFrac1 = mf_rock(Ti .- 273.15)
    end
    return ti, Ti, CrFrac1
end


# function ZirconIsotopeDiffusion()

    Runname = "Test"
    !isdir("Results") && mkpath("Results")

    # parameters for simulations
    CbulkZr = 100
    tyear = 3600*24*365
    iplot = 1 # plot results
    n = 500 # number of points in radial mesh. Can be changed by user depends on desired accuracy
    nt = 500
    CZirc = 490000.0 # zirconium concentration in zircon, ppm
    XH2O = 2 # initial water content in melt, required for diffusion coefficient simulations.
    Tsolidus = 400 + 273 # arbitrary solidus temperature for phase diagram used
    Csupsat = 3 # ppm supersaturation to cause nucleation of a new crystal upon cooling
    UCR = 1 # Critical concentration for microZircon nucleation on major minerals
    ZircNuc = 1e-4 # Zircon stable nuclei in cm
    L = 0.1 # 20e-4*(CZirc/CbulkZr)^(1./3.); radius of melt cell
    DGfZr = 0.5 # ratio of diffusion coefficients of Hf to Zr; change for other element of interest

    # Solve for Tsat
    function equation!(F, T)
        F[1] = ZrSaturation(T[1])*mf_rock(T[1]-273.15) - CbulkZr
    end

    result = nlsolve(equation!, [1000.0])
    Tsat = result.zero[1]
    # global  n R A B D F alpha1 beta Xs UCR ZrPl Tsat  CbulkZr MinCore DGfZr S0
    # global Dplag Temp MeltFrac time XH2O Tsolidus Csupsat V Dscale tscale L S W t CZirc CPl ZircNuc Czl Czh

    # parameters for the simulation
    parameters = Dict(
        "Tsat" => Tsat, # Starting at saturation
        "Tend" => 695,  # final temperature, C
        "tfin" => 1500, # final time
        "Cbulk" => CbulkZr,
        "RhoZrM" => 4.7/2.3, # Ratio of zircon to melt density
        "Kmin" => 0.1,  # Parameters for zirconiun partition coefficient in major phase
        "Kmax" => 0.1,
        "Crit" => 30,
        "delta" => 0.2,
        "Ktrace" => 0.1, # trace partition coefficient in major phase.
        "Trace" => "Hf",
        "XH20" => XH2O,
        "L" => L,
        "DGfZr" => DGfZr,  # diffusion coefficient ratio
        "nt" => nt
    )
    tr = []
    Tr = []
    time, Temp, MeltFrac = TemperatureHistory_m_erupt(tr, Tr, parameters)

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

    # Scaling
    tfin = time[end] # total time in years of the process
    # SCALING-----------------
    Ds = DiffusionCoefficient((750 + 273.15), XH2O, DGfZr)
    Dscale = Ds[1]
    tscale = L^2 / Dscale # dimensionless scale for the time
    time = time ./ tscale
    # nt = L(time)  # this is obsolete as the nt does not change with scaling
    # END:SCALING-----------------

    # Initial Conditions
    t = time[1] / tscale
    ZirconRadius = 2e-4
    Xs = ZirconRadius / L
    ZircNuc = ZircNuc / L
    S = (Xs^3 + MeltFrac[1] * (1 - Xs^3))^(1/3)
    S0 = S
    dt = time[2] - time[1]
    W = 0
    V = 0


    C0[:, 1] .= ZrSaturation(Temp[1]) * 0.5145
    C0[:, 2] .= ZrSaturation(Temp[1]) * 0.1122
    C0[:, 3] .= ZrSaturation(Temp[1]) * 0.1715
    C0[:, 4] .= ZrSaturation(Temp[1]) * 0.1738
    C0[:, 5] .= ZrSaturation(Temp[1]) * 0.0280
    C0[:, 6] .= CZirc / kdHf(Temp[1], parameters) / 70
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
    tt[1] = time[1] * tscale
    Zcomp[1,1] = C0[1, 4] / C0[1, 1]
    ZrHF[1,1] = CZirc / kdHf(Temp[1], parameters) / C0[1, 6]
    # Zcomp[1] = C0[1, 4] / C0[1, 1] #matlab version
    # ZrHF[1] = CZirc / kdHf(Temp[1], par) / C0[1, 6] #matlab version
    Melndelta[1,1] = Zcomp[1,1]

    R = range(0, stop=1, length=n)
    rr = range(0, stop=1, length=n)
    CZircon[1:5] = 4 * π * CZirc * C0[1, 1:5] / ZrSaturation(Temp[1]) * ZirconRadius^3 / 3
    Cplag[1:5] .= 0
    CintS[1, 1:5] = CZircon[1:5] + 4 * π * C0[1, 1:5] * (S^3 - ZirconRadius^3) / 3

    global  n, R,Dplag, ZrPl, MinCore, time, tscale, S0
    global A, B, D, F, alpha, beta, Xs, Temp, MeltFrac, XH2O, Tsolidus, V, W, Csupsat, Dscale, UCR, CZirc, S, ZircNuc, Czl, Czh, Dflux
    # MAIn LOOP in time _______________________

    # Main loop
    for i = 2:nt-1
    # for i = 2:100
        if MeltFrac[i] > 0.01
            C, Czl, Czh, Csat, Dif, S = progonka(C0, dt, i, parameters)
            dt = time[i] - time[i-1]
            C0 = C
        else
            V = 0
            W = 0
        end
        rr = R * (S - Xs) .+ Xs
        Csat = ZrSaturation(Temp[i])
        CZircon[1:5] = CZircon[1:5] - CZirc * C[1, 1:5] / Csat * 4 * π * Xs^2 * V * dt
        Cplag[1:5] = Cplag[1:5] - C[end, 1:5] .* Dplag[1:5] * 4 * π * S^2 * W * dt
        Cint[1:5] .= 0
        for ik = 2:n
            Cint[1:5] = Cint[1:5] + (C[ik-1, 1:5] * rr[ik-1]^2 + C[ik, 1:5] * rr[ik]^2) / 2 * (rr[ik] - rr[ik-1])
        end
        Cint = 4 * π * Cint + parameters["RhoZrM"] * CZircon + Cplag
        CintS[i, 1:5] = Cint[1:5]

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
        tt[i] = time[i] * tscale
        UU[i] = C[1] - ZrSaturation(Temp[i])
        Tsave[i] = Temp[i] - 273
        ZrPls[i] = minimum(XXs[1:i, 1])
        Zcompl[i] = Czl / CZirc
        Zcomph[i] = Czh / CZirc
        Zcomp[1,i] = C[1, 4] / C[1, 1]
        Melndelta[1,i] = Ch / Cl
        ZrHF[i] = CZirc / kdHf(Temp[i], parameters) / C0[1, 6]
        CC[i, 1:n, 1:6] = C0[1:n, 1:6]
    end
    # Plot results (if iplot is set)


    if iplot == 1
        fig = Figure(size = (800, 800), backgroundcolor = :white)

        # Subplots
        ax1 = Axis(fig[1, 1], xlabel = "Time (years)", ylabel = "Zr radius")
        ax2 = Axis(fig[1, 2], xlabel = "Time (years)", ylabel = "Temperature T, C")
        ax3 = Axis(fig[2, 1], xlabel = "Distance", ylabel = L"Growth Rate, μm/a^{-1}")
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
    # Save results

    # Print the figure to a PDF file
    CairoMakie.save("Results/Test.pdf", fig)

    i = nt - 2

    # Convert array to DataFrame
    Rsave = DataFrame(time_ka = tt[1:i-1] / 1e3, Rad_um = XXs[1:i-1], Gr_rate_mm_a = VV[1:i-1], Temp_C = Tsave[1:i-1], DelZr = DelZr[1:i-1], DelMlt = DelMlt[1:i-1], ZrHf = ZrHF[1:i-1])

    # Write DataFrame to CSV file
    CSV.write("Results/$Runname.csv", Rsave)

    # Append structure to CSV file
    par = DataFrame(fname = "$Runname")

    if !isfile("Results/summary.csv")
        wwar = true
    else
        wwar = false
    end

    if wwar
        CSV.write("Results/summary.csv", par)
    else
        existing = CSV.read("Results/summary.csv", DataFrame)
        append!(existing, par)
        CSV.write("Results/summary.csv", existing)
    end
# end

#Matlab
# >> sum(C(1,:))
# ans = 129.4123
@test sum(C[1,:]) ≈ 129.4123 atol=1e-4

# >> max(XXs)
# ans = 19.4589
@test maximum(XXs) ≈ 19.4589 atol=1e-4

# >> min(DelZr)
# ans = -2.0528
# Julia = -2.5351 (tolerance??)
@test minimum(DelZr[2:end-1]) ≈ -2.0528 atol=1e-3
