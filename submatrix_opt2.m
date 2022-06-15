function [shiftSpec, shiftConc, CalcPara_submatrix] = submatrix_opt2(matrix, CalcPara_submatrix);

%     sig_size = CalcPara_submatrix.sig_size;
%     x_range = CalcPara_submatrix.x_range;
%     x = linspace(x_range(1), x_range(2), sig_size);
%     calc_range = CalcPara_submatrix.calc_range;
%     calc_idx = find(x >= calc_range(1) & x <= calc_range(2));
    nit = CalcPara_submatrix.NumIter;
    calc_idx = CalcPara_submatrix.calc_idx;
    matrix = matrix(:, calc_idx);
%     matrix = matrix';
    % addition of noise comp
%     n_comp = CalcPara_submatrix.n_comp + 1;
    n_comp = CalcPara_submatrix.n_comp;
    [conc, inispec] = nnmf(matrix, n_comp);

    [szcomp, szsig] = size(inispec);

    %**************************************************************************
    % MCR ALS OPTIMIZATION,
    %**************************************************************************

    % A) DATA PREPARATION AND INPUT

    % INITIALIZATIONS
    SpecNonneg = CalcPara_submatrix.Spectral_Nonnegativity; % on / off  ...If on, apply nonnegative constratint
    ConcNonneg = CalcPara_submatrix.Concentration_Nonnegativity; % on / off  ....If on,  apply nonnegative constratint
%     SpecNonposi = CalcPara_submatrix.Spectral_Nonpositivity; %ex [0;1;0]
    alsOptions = CalcPara_submatrix.alsOptions;
    
    nc = min(size(inispec));
    % non-negativity setting
    if ConcNonneg == "on"
        alsOptions.nonegC.noneg = 1;
        alsOptions.nonegC.cneg = ones(1, nc);
    elseif ConcNonneg == "off"
        alsOptions.nonegC.noneg = 0;
        alsOptions.nonegC.cneg = zeros(1, nc);
    end

    if SpecNonneg == "on"
        alsOptions.nonegS.noneg = 1;
        alsOptions.nonegS.spneg = ones(nc, 1);
    elseif SpecNonneg == "off"
        alsOptions.nonegC.noneg = 0;
        alsOptions.nonegS.spneg = zeros(nc, 1);
    end
    
    alsOptions.resultats.optim_niter = nit;
    alsOptions.opt.nit = nit;
    
    matdad = matrix;
    iniesta = inispec;
    nsign = min(size(iniesta));
    [nrow, ncol] = size(matdad);
    [nrow2, ncol2] = size(iniesta);
    ils = 2; % need to add
    nexp = 1;
%     isp = ones(nsign, 1);

    if ils == 2,
        abss = iniesta;
        [nsign, ncol] = size(abss);
    end

    if nexp == 1,
        ncinic(nexp) = 1;
        ncfin(nexp) = ncol;
        nrinic(nexp) = 1;
        nrfin(nexp) = nrow;
    end

    niter = 0; % iterations counter
    idev = 0; % divergence counter
    idevmax = 10; % maximum number of diverging iterations
    answer = 'n'; % default answer
    ineg = 0; % used for non-negativity constraints
    imod = 0; % used for unimodality constraints
    iclos0 = 0; % used for closure constraints
    datamod = 99; % in three-way type of matrix augmentation (1=row, 2=column)
    vclos1n = 0; % used for closure constraints
    vclos2n = 0; % used for closure constraints
    inorm = 0; % no normalizatio (when closurfe is applied)

    %***************************
    % DEFINITION OF THE DATA SET
    %***************************

    totalconc(1:nsign, 1:nexp) = ones(nsign, nexp);

    % WHEN ONLY ONE EXPERIMENT IS PRESENT EVERYTHING IS CLEAR

    if nexp == 1
        nrsol(1) = nrow;
        nrinic(1) = 1;
        nrfin(1) = nrsol(1);
        nesp(1) = nsign;
        matr = 1;
        matc = 1;
        isp(1, 1:nsign) = ones(1, nsign);
    end
    ishape = 0;

    % ***********************************************************
    % B) REPRODUCTION OF THE ORIGINAL DATA MATRIX BY PCA
    % ***********************************************************

    % dn is the experimental matrix and d is the PCA reproduced matrix
    d = matdad;
    dn = d;
%     [u, s, v, d, sd] = pcarep(dn, nsign);
    [u, s, v] = svd(dn, 0);
    u = u(:, 1:nsign);
    s = s(1:nsign, 1:nsign);
    v = v(:, 1:nsign);
    d = u * s * v';

    sstn = sum(sum(dn .* dn));
    sst = sum(sum(d .* d));
    sigma2 = sqrt(sstn);

    %save('temp.mat');
    % **************************************************************
    % C) STARTING ALTERNATING CONSTRAINED LEAST SQUARES OPTIMIZATION
    % **************************************************************
    
    plot_lof = [];
    plot_R2 = [];
    plot_sigmaC = [];
    plot_sigmaN = [];
    plot_specs = [];
    plot_concs = [];
    
    while niter < alsOptions.opt.nit;

        niter = niter + 1;

        % ******************************************
        % D) ESTIMATE CONCENTRATIONS (ALS solutions)
        % ******************************************

        conc = d / abss;

        % ******************************************
        % CONSTRAIN APPROPRIATELY THE CONCENTRATIONS
        % ******************************************

        % ****************
        % non-negativity
        % ****************

        if alsOptions.nonegC.noneg == 1;
            ineg = alsOptions.nonegC.noneg;
            if ineg == 1;

                ialg = alsOptions.nonegC.ialg;
                ncneg = alsOptions.nonegC.ncneg;
                cneg = alsOptions.nonegC.cneg;

                for i = 1:matc
                    kinic = nrinic(i);
                    kfin = nrfin(i);
                    conc2  = conc(kinic:kfin, :);

                    if ialg == 0
                        for k = 1:nsign,
                            if cneg(i, k) == 1
                                for j = 1:kfin+1-kinic,
                                    if conc2(j,k) < 0.0,
                                        conc2(j,k) = 0.0;
                                    end
                                end
                            end
                        end
                    end

                    if ialg == 1
                        for j = kinic:kfin
                            if cneg(i, :) == ones(1, size(isp, 2))
                                x = lsqnonneg(abss', d(j,:)');
                                conc2(j - kinic + 1, :) = x';
                            end
                        end
                    end

                    if ialg == 2
                        for j = kinic:kfin
                            if cneg(i,:) == ones(1, size(isp, 2))
                                x = fnnls(abss * abss', abss * d(j,:)');
                                conc2(j - kinic + 1,:) = x';
                            end
                        end
                    end

                    conc(kinic:kfin, :) = conc2;
                end
            end
        end

        ishape=0;

        if ishape == 0 | niter == 1,
            for j = 1:nsign,
                for inexp = 1:matc,
                    totalconc(j, inexp) = sum(conc(nrinic(inexp):nrfin(inexp), j));
                end
                if totalconc(j, 1) > 0,
                    rt(j, 1:matc) = totalconc(j, 1:matc) ./ totalconc(j, 1);
                else
                    rt(j, 1:matc) = totalconc(j, 1:matc);
                end
            end
        end

        % areas under concentration profiles
        area = totalconc;

        % ********************************
        % ESTIMATE SPECTRA (ALS solution)
        % ********************************

        abss=conc \ d;


        % ********************
        % non-negative spectra
        % ********************

        ineg = alsOptions.nonegS.noneg;
        if ineg == 1,

            ialgs = alsOptions.nonegS.ialgs;
            nspneg = alsOptions.nonegS.nspneg;
            spneg = alsOptions.nonegS.spneg;
            if matr > 1
                ncinic = alsOptions.multi.ncinic;
                ncfin  = alsOptions.multi.ncfin;
            end

            for i = 1:matr
                kinic = ncinic(i);
                kfin = ncfin(i);
                abss2 = abss(:, kinic:kfin);

                if ialgs == 0,
                    for k = 1:nsign,
                        if spneg(k, i) == 1
                            for j = 1:kfin+1-kinic,
                                if abss2(k, j) < 0.0,
                                    abss2(k, j) = 0.0;
                                end
                            end
                        end
                    end
                end

                if ialgs == 1,
                    for j = kinic:kfin,
                        if spneg(:, i) == ones(size(isp, 2), 1)
                            abss2(:, j - kinic + 1) = lsqnonneg(conc, d(:, j));
                        end
                    end
                end

                if ialgs == 2,
                    for j = kinic:kfin,
                        if spneg(:, i) == ones(size(isp, 2), 1)
                            abss2(:, j - kinic + 1) = fnnls(conc' * conc, conc' * d(:,j));
                        end
                    end
                end
                abss(:, kinic:kfin) = abss2;
            end
        end    

        % invert sign of some components 
%         if sum(SpecNonposi) ~= 0
%             idx_negative = SpecNonposi == 1;
%             abss(idx_negative, :) = - abss(idx_negative, :);
%         end
    %         
        % ************************
        % NORMALIZATION OF SPECTRA
        % ************************

        % equal height
        inorm = alsOptions.closure.inorm;

        if inorm == 1;
            maxabss = max(abss');
            for i = 1:nsign,
                abss(i, :) = abss(i, :) ./ maxabss(i);
            end
        end

        % equal length - divided by Frobenius Norm

        if inorm == 2, abss  = normv2(abss); end

        % equal length - divided by Total Sum Norm

        if inorm == 3, abss = normv3(abss); end           
    
        % *******************************
        % CALCULATE RESIDUALS
        % *******************************

        res= d - conc * abss;
        resn = dn - conc * abss;

        % ********************************
        % OPTIMIZATION RESULTS
        % *********************************

        disp(' ' ); disp(' '); disp(['ITERATION ', num2str(niter)]);
        u = sum(sum(res .* res));
        un = sum(sum(resn .* resn));
        disp(['Sum of squares respect PCA reprod. = ', num2str(u)]);
        sigma = sqrt(u / (nrow * ncol));
        sigman  =sqrt(un / (nrow * ncol));
        disp(['Old sigma = ', num2str(sigma2),' -----> New sigma = ', num2str(sigma)]);
        disp(['Sigma respect experimental data = ', num2str(sigman)]);
        disp(' ');
        change =((sigma2 - sigma) / sigma);  % sigma2 is previous sigma

        if change < 0.0,
            disp(' ')
            disp('FITING IS NOT IMPROVING !!!')
            idev = idev + 1;
        else,
            disp('FITING IS IMPROVING !!!')
            idev = 0;
        end

        change = change * 100;
        disp(['Change in sigma (%) = ', num2str(change)]);
        sstd(1) = sqrt(u / sst) * 100;  % sst and sstn are normalization constant
        sstd(2) = sqrt(un / sstn) * 100;
        disp(['Fitting error (lack of fit, lof) in % (PCA) = ', num2str(sstd(1))]);
        disp(['Fitting error (lack of fit, lof) in % (exp) = ', num2str(sstd(2))]);
        r2 = (sstn - un) / sstn;
        disp(['Percent of variance explained (r2) is ', num2str(100 * r2)]);

        % save parameters
        % ******************************************** used for plotting
        % optimization information
        % *****************************************************************

        plot_lof=[plot_lof; sstd(2)];
        plot_sigmaC = [plot_sigmaC; change];
        plot_sigmaN = [plot_sigmaN; sigman];
        plot_R2 = [plot_R2; 100 * r2];
        plot_specs = [plot_specs; abss];
        plot_concs = [plot_concs conc];


        % ********************
        % DISPLAY PURE SPECTRA
        % ********************
        
        figure(2)
        subplot(2, 1, 1)
        plot(abss')
        drawnow
        subplot(2, 1, 2)
        plot(conc)
        drawnow

%             clear cx_plot sx_plot niter_plot change_plot sstd_plot als_end
%         end

        % *************************************************************
        % If change is positive, the optimization is working correctly
        % *************************************************************

        if change > 0 | niter == 1,

            sigma2 = sigma;
            copt = conc;
            sopt = abss;
            sdopt = sstd;
            ropt = res;
            rtopt = rt';
            itopt = niter;
            areaopt = area;
            r2opt = r2;
        end

        % ******************************************************************
        % test for convergence within maximum number of iterations allowed
        % ******************************************************************

        if abs(change) <  alsOptions.opt.tolsigma,

            % finish the iterative optimization because convergence is achieved

            disp(' '); disp(' ');
            disp('CONVERGENCE IS ACHIEVED !!!!')
            disp(' ')
            disp(['Fitting error (lack of fit, lof) in % at the optimum = ', num2str(sdopt(1, 1)), '(PCA) ', num2str(sdopt(1, 2)), '(exp)']);
            disp(['Percent of variance explained (r2)at the optimum is ', num2str(100 * r2opt)]);
            disp('Relative species conc. areas respect matrix (sample) 1at the optimum'), disp(rtopt')
            disp(['Plots are at optimum in the iteration ', num2str(itopt)]);


            if isempty(alsOptions.out.out_conc) == 0
                alsOptions.out.out_conc = copt;
            end
            if isempty(alsOptions.out.out_spec) == 0
                alsOptions.out.out_spec = sopt;
            end
            if isempty(alsOptions.out.out_res) == 0
                alsOptions.out.out_res = ropt;
            end
            if isempty(alsOptions.out.out_std) == 0
                alsOptions.out.out_std = sdopt;
            end
            if isempty(alsOptions.out.out_area) == 0
                alsOptions.out.out_area = areaopt;
            end
            if isempty(alsOptions.out.out_rat) == 0
                alsOptions.out.out_rat = rtopt;
            end
%             return         % 1st return (end of the optimization, convergence)
        end

        %  finish the iterative optimization if divergence occurs 20 times consecutively
        if idev > 20,
            disp(' ');disp(' ');
            disp('FIT NOT IMPROVING FOR 20 TMES CONSECUTIVELY (DIVERGENCE?), STOP!!!')
            disp(' ')
            disp(['Fitting error (lack of fit, lof) in % at the optimum = ', num2str(sdopt(1, 1)),'(PCA) ', num2str(sdopt(1, 2)), '(exp)']);
            disp(['Percent of variance explained (r2)at the optimum is ', num2str(100 * r2opt)]);
            disp('Relative species conc. areas respect matrix (sample) 1 at the optimum'), disp(rtopt)
            disp(['Plots are at optimum in the iteration ', num2str(itopt)]);


            if isempty(alsOptions.out.out_conc) == 0
                alsOptions.out.out_conc = copt;
            end
            if isempty(alsOptions.out.out_spec) == 0
                alsOptions.out.out_spec = sopt;
            end
            if isempty(alsOptions.out.out_res) == 0
                alsOptions.out.out_res = ropt;
            end
            if isempty(alsOptions.out.out_std) == 0
                alsOptions.out.out_std = sdopt;
            end
            if isempty(alsOptions.out.out_area) == 0
                alsOptions.out.out_area = areaopt;
            end
            if isempty(alsOptions.out.out_rat) == 0
                alsOptions.out.out_rat = rtopt;
            end

%             return          % 2nd return (end of optimization, divergence)        
        end

        % this end refers to number of iterations initially proposed exceeded
        if niter == alsOptions.opt.nit
            % finish the iterative optimization if maximum number of allowed iterations is exceeded
            disp(' ');disp(' ');
            disp('NUMBER OF ITERATIONS EXCEEDED THE ALLOWED!')
            disp(' ')
            disp(['Fitting error (lack of fit, lof) in % at the optimum = ', num2str(sdopt(1, 1)),'(PCA) ', num2str(sdopt(1, 2)), '(exp)']);
            disp(['Percent of variance explained (r2)at the optimum is ',num2str(100 * r2opt)]);
            disp('Relative species conc. areas respect matrix (sample) 1 at the optimum'), disp(rtopt)
            disp(['Plots are at optimum in the iteration ', num2str(itopt)]);

            if isempty(alsOptions.out.out_conc) == 0
                alsOptions.out.out_conc = copt;
            end
            if isempty(alsOptions.out.out_spec) == 0
                alsOptions.out.out_spec = sopt;
            end
            if isempty(alsOptions.out.out_res) == 0
                alsOptions.out.out_res = ropt;
            end
            if isempty(alsOptions.out.out_std) == 0
                alsOptions.out.out_std = sdopt;
            end
            if isempty(alsOptions.out.out_area) == 0
                alsOptions.out.out_area = areaopt;
            end
            if isempty(alsOptions.out.out_rat) == 0
                alsOptions.out.out_rat = rtopt;
            end

%             return      %3rd return (end of optimization, number of iterations exceeded)
        end
    end
    
    alsOptions.resultats.plot_lof = plot_lof;
    alsOptions.resultats.plot_R2 = plot_R2;
    alsOptions.resultats.plot_sigmaC = plot_sigmaC;
    alsOptions.resultats.plot_sigmaN = plot_sigmaN;
    alsOptions.resultats.total_niter = niter;
    alsOptions.resultats.plot_specs = plot_specs;
    alsOptions.resultats.plot_concs = plot_concs;
    
    shiftSpec = sopt;
    shiftConc = copt;
    % sub_conc = mean(copt,2);
    % sub_abss = mean(sopt);
    % do not divide both, it will double the effect.
%     sub_conc = sum(copt, 2);
%     sub_abss = mean(sopt);
    
    CalcPara_submatrix.alsOptions = alsOptions;

end
