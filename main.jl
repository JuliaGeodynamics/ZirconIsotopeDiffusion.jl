using NLsolve
using CairoMakie
using Interpolations
using Trapz
using CSV
using DataFrames

function ZirconIsotopeDiffusion()

    include("src/Helperfunctions.jl")
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
    length = 0.1 # 20e-4*(CZirc/CbulkZr)^(1./3.); radius of melt cell
    DGfZr = 0.5 # ratio of diffusion coefficients of Hf to Zr; change for other element of interest

    # Solve for Tsat
    function equation!(F, T)
        F[1] = ZrSaturation(T[1])*mf_rock(T[1]-273.15) - CbulkZr
    end

    result = nlsolve(equation!, [1000.0])
    Tsat = result.zero[1]
    # global  n R A B D F alpha1 beta Xs UCR ZrPl Tsat  CbulkZr MinCore DGfZr S0
    # global Dplag Temp MeltFrac time XH2O Tsolidus Csupsat V Dscale tscale length S W t CZirc CPl ZircNuc Czl Czh

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
        "L" => length,
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
    tscale = length^2 / Dscale # dimensionless scale for the time
    time = time ./ tscale
    # nt = length(time)  # this is obsolete as the nt does not change with scaling
    # END:SCALING-----------------

    # Initial Conditions
    t = time[1] / tscale
    ZirconRadius = 2e-4
    Xs = ZirconRadius / length
    ZircNuc = ZircNuc / length
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
    XXs[1] = Xs * 1e4 * length
    RRd[1] = S * 1e4 * length
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
            lines!(ax1, rr * length * 1e4, (C[:, 4] .* 0.5145 ./ C[:, 1] ./ 0.1738 .- 1) .* 1000, linewidth = 1.5)
            # lines!(ax2, rr * length * 1e4, sum(C[:, 1:5], dims = 2) ./ C[:, 6], linewidth = 1.5)
            display(fig)
        end


        t += dt
        rr = R * (S - Xs) .+ Xs
        Cl = trapz(rr[:,1], rr[:,1].^2 .* C[:, 1])
        Ch = trapz(rr, rr.^2 .* C[:, 4])
        Xs = max(ZircNuc, Xs - V * dt)
        S0 = S
        XXs[i] = Xs * 1e4 * length  # zircon radius in um
        RRd[i] = S * 1e4 * length  # melt cell radius in um
        VV[i] = -V * length * 1e4 / tscale  # array of dissolution rate
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
end

ZirconIsotopeDiffusion()
