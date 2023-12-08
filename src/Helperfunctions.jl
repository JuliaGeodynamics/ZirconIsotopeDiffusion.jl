#Helperfunctions


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

function mf_rock(T)
    t2 = T .* T
    t7 = exp.(0.961026371384066e3 .- 0.3590508961e1 .* T .+ 0.4479483398e-2 .* t2 .- 0.1866187556e-5 .* t2 .* T)
    CF = 0.1e1 ./ (0.1e1 .+ t7)
    return CF
end

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
