function [sub_conc, sub_abss, shiftSig, shiftConc] = submatrix_opt(matrix, CalcPara_submatrix);

    sig_size = CalcPara_submatrix.sig_size;
    x_range = CalcPara_submatrix.x_range;
    x = linspace(x_range(1), x_range(2), sig_size);
    calc_range = CalcPara_submatrix.calc_range;
    calc_idx = find(x >= calc_range(1) & x <= calc_range(2));
    matrix = matrix(calc_idx, :);
    matrix = matrix';
    n_comp = 3;
    [conc, inispec] = nnmf(matrix, n_comp);

    [szcomp, szsig] = size(inispec);

%**************************************************************************
% MCR ALS OPTIMIZATION,
%**************************************************************************

% A) DATA PREPARATION AND INPUT

% INITIALIZATIONS
SpecNonneg = CalcPara_submatrix.Spectral_Nonnegativity; % on / off  ...If on, apply nonnegative constratint
ConcNonneg = CalcPara_submatrix.Concentration_Nonnegativity; % on / off  ....If on,  apply nonnegative constratint
SpecNonposi = CalcPara_submatrix.Spectral_Nonpositivity; %ex [0;1;0]
mcr_als = CalcPara_submatrix.mcr_als;
matdad = matrix;
iniesta = inispec;
nsign = min(size(iniesta));
[nrow, ncol] = size(matdad);
[nrow2, ncol2] = size(iniesta);
ils = 2; % need to add
nexp = 1;
    isp = ones(nsign, 1);

if nrow2==nrow,	nsign=ncol2;
    ils=1;
end
if ncol2==nrow, nsign=nrow2; 
    iniesta=iniesta'; 
    ils=1; 
end

if ncol2==ncol, nsign=nrow2; 
    ils=2;
end
if nrow2==ncol, nsign=ncol2; 
    iniesta=iniesta'; 
    ils=2; 
end

if ils==1,
    conc=iniesta;
    [nrow,nsign]=size(conc);
    abss=conc\matdad;
end

if ils==2,
    abss=iniesta;
    [nsign,ncol]=size(abss);
end

if nexp==1,
    ncinic(nexp)=1;
    ncfin(nexp)=ncol;
    nrinic(nexp)=1;
    nrfin(nexp)=nrow;
end

niter=0;% iterations counter
idev=0;% divergence counter
idevmax=10;% maximum number of diverging iterations
answer='n'; % default answer
ineg=0;% used for non-negativity constraints
imod=0;% used for unimodality constraints
iclos0=0;% used for closure constraints
iassim=0;% used for shape constraints
datamod=99;% in three-way type of matrix augmentation (1=row,2=column)
vclos1n=0;% used for closure constraints
vclos2n=0;% used for closure constraints
inorm=0;% no normalizatio (when closurfe is applied)
type_csel=[]; %no equality/lower than constraints in concentrations
type_ssel=[]; %no equality/lower than constraints in spectra

%***************************
% DEFINITION OF THE DATA SET
%***************************

totalconc(1:nsign,1:nexp)=ones(nsign,nexp);

% WHEN ONLY ONE EXPERIMENT IS PRESENT EVERYTHING IS CLEAR

if nexp==1
    nrsol(1)=nrow;
    nrinic(1)=1;
    nrfin(1)=nrsol(1);
    nesp(1)=nsign;
    matr = 1;
    matc = 1;
    isp(1,1:nsign)=ones(1,nsign);
end
ishape=0;

% ***********************************************************
% B) REPRODUCTION OF THE ORIGINAL DATA MATRIX BY PCA
% ***********************************************************

% dn is the experimental matrix and d is the PCA reproduced matrix
d = matrix;
dn=d;
[u,s,v,d,sd]=pcarep(dn,nsign);

sstn=sum(sum(dn.*dn));
sst=sum(sum(d.*d));
sigma2=sqrt(sstn);

%save('temp.mat');
% **************************************************************
% C) STARTING ALTERNATING CONSTRAINED LEAST SQUARES OPTIMIZATION
% **************************************************************

while niter < mcr_als.alsOptions.opt.nit ;
    
    niter=niter+1;
    
    % ******************************************
    % D) ESTIMATE CONCENTRATIONS (ALS solutions)
    % ******************************************
    
    conc=d/abss;
    
    % ******************************************
    % CONSTRAIN APPROPRIATELY THE CONCENTRATIONS
    % ******************************************
    
    % ****************
    % non-negativity
    % ****************
    
    if mcr_als.alsOptions.nonegC.noneg==1;
        ineg=mcr_als.alsOptions.nonegC.noneg;
        if ineg ==1;
            
            ialg=mcr_als.alsOptions.nonegC.ialg;
            ncneg=mcr_als.alsOptions.nonegC.ncneg;
            cneg=mcr_als.alsOptions.nonegC.cneg;
            
            for i=1:matc
                kinic=nrinic(i);
                kfin=nrfin(i);
                conc2=conc(kinic:kfin,:);
                
                if ialg==0
                    for k=1:nsign,
                        if cneg(i,k) ==1
                            for j=1:kfin+1-kinic,
                                if conc2(j,k)<0.0,
                                    conc2(j,k)=0.0;
                                end
                            end
                        end
                    end
                end
                
                if ialg==1
                    for j=kinic:kfin
                        if cneg(i,:) == ones(1,size(isp,2))
                            x=lsqnonneg(abss',d(j,:)');
                            conc2(j-kinic+1,:)=x';
                        end
                    end
                end
                
                if ialg==2
                    for j=kinic:kfin
                        if cneg(i,:) == ones(1,size(isp,2))
                            x=fnnls(abss*abss',abss*d(j,:)');
                            conc2(j-kinic+1,:)=x';
                        end
                    end
                end
                
                conc(kinic:kfin,:) = conc2;
            end
        end
    end

    % **************************
    % zero concentration species
    % **************************
    
    if matc>1
        for i=1:matc,
            for j=1:nsign,
                if isp(i,j)==0,
                    conc(nrinic(i):nrfin(i),j)=zeros(nrsol(i),1);
                end
            end
        end
    end
    
        
    % ************************************************
    % QUANTITATIVE INFORMATION FOR THREE-WAY DATA SETS
    % ************************************************
    
    if evalin('base','mcr_als.alsOptions.trilin.appTril')==1
        ishape=evalin('base','mcr_als.alsOptions.trilin.ishape');
    else
        ishape=0;
    end
    
    if ishape==0 | niter==1,
        for j=1:nsign,
            for inexp=1:matc,
                totalconc(j,inexp)=sum(conc(nrinic(inexp):nrfin(inexp),j));
            end
            if totalconc(j,1)>0,
                rt(j,1:matc)=totalconc(j,1:matc)./totalconc(j,1);
            else
                rt(j,1:matc)=totalconc(j,1:matc);
            end
        end
    end
    
    % areas under concentration profiles
    area=totalconc;
    
    % ********************************
    % ESTIMATE SPECTRA (ALS solution)
    % ********************************
    
    abss=conc\d;
    
    
    % ********************
    % non-negative spectra
    % ********************
    
    ineg=mcr_als.alsOptions.nonegS.noneg;
    if ineg==1,
        
        ialgs=mcr_als.alsOptions.nonegS.ialgs;
        nspneg=mcr_als.alsOptions.nonegS.nspneg;
        spneg=mcr_als.alsOptions.nonegS.spneg;
        if matr>1
            ncinic=mcr_als.alsOptions.multi.ncinic;
            ncfin=mcr_als.alsOptions.multi.ncfin;
        end
        
        for i = 1:matr
            kinic = ncinic(i);
            kfin = ncfin(i);
            abss2 = abss(:,kinic:kfin);
            
            if ialgs==0,
                for k=1:nsign,
                    if spneg(k,i)==1
                        for j=1:kfin+1-kinic,
                            if abss2(k,j)<0.0,
                                abss2(k,j)=0.0;
                            end
                        end
                    end
                end
            end
            
            if ialgs==1,
                for j=kinic:kfin,
                    if spneg(:,i)== ones(size(isp,2),1)
                        abss2(:,j-kinic+1)=lsqnonneg(conc,d(:,j));
                    end
                end
            end
            
            if ialgs==2,
                for j=kinic:kfin,
                    if spneg(:,i)== ones(size(isp,2),1)
                        abss2(:,j-kinic+1)=fnnls(conc'*conc,conc'*d(:,j));
                    end
                end
            end
            abss(:,kinic:kfin)=abss2;
        end
    end
    
    
    %     % invert sign of some components 
    if sum(SpecNonposi) ~= 0
        idx_negative = SpecNonposi == 1;
        abss(idx_negative,:) = - abss(idx_negative,:);
    end
%     
    
%     %%%%%%%%%%%%%
%     %Weight constratints 
%     
%     if isempty(CalcPara_submatrix.Opt.Weight)==0
%         abss_original = iniesta;
%         
%         if idev > NumIterDivergence
%             clear
%             idstrgrest='on';
%             load('temp.mat')
%             continue
%         end
%         
%         if exist('idstrgrest')==1
%             
%             freeline = mcr_weight.freeline;
%             allowance = mcr_weight.allowance;
%             size_abss = size(abss);
%             if mcr_weight.limit == "on"
%                 
%                 abss =weight.*abss;
%                 area=sum(abss,2);
%                 area_ini=sum(inispec,2);
%                 
%                 if any(area < Area_Constraints(1) * area_ini) || any(area > Area_Constraints(2) * area_ini)
%                     Locs_A=find(area < Area_Constraints(1) * area_ini | area > Area_Constraints(2) * area_ini );
%                     abss(Locs_A,:)=inispec(Locs_A,:);
%                 end
%             end
%             
%         else
%             abss=weight.*abss;
%             
%         end
%         
%     end
    
    
    % ************************
    % NORMALIZATION OF SPECTRA
    % ************************
    
    
    if evalin('base','mcr_als.alsOptions.closure.closure')==0
        
        % equal height
        if evalin('base','mcr_als.alsOptions.correlation.checkSNorm')==0;
            inorm=evalin('base','mcr_als.alsOptions.closure.inorm');
            
            if inorm==1;
                maxabss=max(abss');
                for i=1:nsign,
                    abss(i,:)=abss(i,:)./maxabss(i);
                end
            end
            
            % equal length - divided by Frobenius Norm
            
            if inorm==2, abss=normv2(abss); end
            
            % equal length - divided by Total Sum Norm
            
            if inorm==3, abss=normv3(abss); end
            
            % in case of application of a correlation constraint,
            % spectra of non-correlated components can go out of scale
            % and produce rank deficiency situations
            % They should be in the same scale as the one with the correlation constraint
        elseif evalin('base','mcr_als.alsOptions.correlation.checkSNorm')==1;
            for i=1:nsign,
                if compreg(i)==1,
                    maxabss=max(abss(i,:)');
                    imax=i;
                end
            end
            for i=1:nsign,
                % if compreg(i)==0 | i~=imax,
                if compreg(i)==0,
                    scaleabss=maxabss/max(abss(i,:)');
                    abss(i,:)=abss(i,:).*scaleabss;
                end
            end
            
        end        
    end
end
    
   
    % *******************************
    % CALCULATE RESIDUALS
    % *******************************
    
    res=d-conc*abss;
    resn=dn-conc*abss;
    
    % ********************************
    % OPTIMIZATION RESULTS
    % *********************************
    
    disp(' ' );disp(' ');disp(['ITERATION ',num2str(niter)]);
    u=sum(sum(res.*res));
    un=sum(sum(resn.*resn));
    disp(['Sum of squares respect PCA reprod. = ', num2str(u)]);
    sigma=sqrt(u/(nrow*ncol));
    sigman=sqrt(un/(nrow*ncol));
    disp(['Old sigma = ', num2str(sigma2),' -----> New sigma = ', num2str(sigma)]);
    disp(['Sigma respect experimental data = ', num2str(sigman)]);
    disp(' ');
    change=((sigma2-sigma)/sigma);
    
    if change < 0.0,
        disp(' ')
        disp('FITING IS NOT IMPROVING !!!')
        idev=idev+1;
    else,
        disp('FITING IS IMPROVING !!!')
        idev=0;
    end
    
    change=change*100;
    disp(['Change in sigma (%) = ', num2str(change)]);
    sstd(1)=sqrt(u/sst)*100;
    sstd(2)=sqrt(un/sstn)*100;
    disp(['Fitting error (lack of fit, lof) in % (PCA) = ', num2str(sstd(1))]);
    disp(['Fitting error (lack of fit, lof) in % (exp) = ', num2str(sstd(2))]);
    r2=(sstn-un)/sstn;
    disp(['Percent of variance explained (r2) is ',num2str(100*r2)]);
    
    % save parameters
    % ******************************************** used for plotting
    % optimization information
    % *****************************************************************
    plot_lof=[];
plot_R2=[];
plot_sigmaC=[];
plot_sigmaN=[];
plot_specs=[];
plot_concs=[];
% kinetic constraint
save_kopt=[];
save_sigK=[];
save_ssq=[];
% correlation constraint
save_yout=[];
save_ycal=[];
save_stats=[];

    plot_lof=[plot_lof;sstd(2)];
    plot_sigmaC=[plot_sigmaC;change];
    plot_sigmaN=[plot_sigmaN;sigman];
    plot_R2=[plot_R2;100*r2];
    plot_specs=[plot_specs;abss];
    plot_concs=[plot_concs conc];
    
    assignin('base','plot_lof',plot_lof);
    assignin('base','plot_R2',plot_R2);
    assignin('base','plot_sigmaC',plot_sigmaC);
    assignin('base','plot_sigmaN',plot_sigmaN);
    assignin('base','total_niter',niter);
    assignin('base','plot_specs',plot_specs);
    assignin('base','plot_concs',plot_concs);
    
    evalin('base','mcr_als.alsOptions.resultats.plot_lof=plot_lof;');
    evalin('base','mcr_als.alsOptions.resultats.plot_R2=plot_R2;');
    evalin('base','mcr_als.alsOptions.resultats.plot_sigmaC=plot_sigmaC;');
    evalin('base','mcr_als.alsOptions.resultats.plot_sigmaN=plot_sigmaN;');
    evalin('base','mcr_als.alsOptions.resultats.total_niter=total_niter;');
    evalin('base','mcr_als.alsOptions.resultats.plot_specs=plot_specs;');
    evalin('base','mcr_als.alsOptions.resultats.plot_concs=plot_concs;');
    
    %     evalin('base','clear plot_lof plot_R2 total_niter plot_specs plot_concs plot_sigmaC plot_sigmaN');
    
    appCinetic=evalin('base','mcr_als.alsOptions.kinetic.appKinetic');
    if appCinetic==1;
        
        kiter=evalin('base','mcr_als.alsOptions.kinetic.results.kopt;');
        sigKiter=evalin('base','mcr_als.alsOptions.kinetic.results.sig_knglm;');
        ssqx=evalin('base','mcr_als.alsOptions.kinetic.results.ssq;');
        
        save_ssq=[save_ssq;ssqx];
        save_kopt=[save_kopt;kiter];
        save_sigK=[save_sigK;sigKiter];
        
        assignin('base','save_kopt',save_kopt);
        assignin('base','save_sigK',save_sigK);
        assignin('base','save_ssq',save_ssq);
        
        evalin('base','mcr_als.alsOptions.kinetic.results.save_kopt=save_kopt;');
        evalin('base','mcr_als.alsOptions.kinetic.results.save_sigK=save_sigK;');
        evalin('base','mcr_als.alsOptions.kinetic.results.save_ssq=save_ssq;');
        
        evalin('base','clear save_kopt save_sigK save_ssq');
        
    end
    
    appCorrel=evalin('base','mcr_als.alsOptions.correlation.appCorrelation;');
    if appCorrel==1;
        
        % yout,ycal, stats available
        
        save_yout=[save_yout;yout];
        save_ycal=[save_ycal;ycal];
        save_stats=[save_stats;stats];
        
        assignin('base','save_yout',save_yout);
        assignin('base','save_ycal',save_ycal);
        assignin('base','save_stats',save_stats);
        
        evalin('base','mcr_als.alsOptions.correlation.results.save_yout=save_yout;');
        evalin('base','mcr_als.alsOptions.correlation.results.save_ycal=save_ycal;');
        evalin('base','mcr_als.alsOptions.correlation.results.save_stats=save_stats;');
        
        evalin('base','clear save_yout save_ycal save_stats');
        
        
    end
    
    
    % ********************
    % DISPLAY PURE SPECTRA
    % ********************
    
    if  evalin('base','mcr_als.alsOptions.opt.gr')=='y',
        
        als_end=0;
        
        assignin('base','cx_plot',conc);
        assignin('base','sx_plot',abss);
        assignin('base','niter_plot',niter);
        assignin('base','change_plot',change);
        assignin('base','sstd_plot',sstd);
        assignin('base','als_end',als_end);
        als_res;
        
        evalin('base','clear cx_plot sx_plot niter_plot change_plot sstd_plot als_end');
    end
    
    
    % *************************************************************
    % If change is positive, the optimization is working correctly
    % *************************************************************
    
    if change>0 | niter==1,
        
        sigma2=sigma;
        copt=conc;
        sopt=abss;
        sdopt=sstd;
        ropt=res;
        rtopt=rt';
        itopt=niter;
        areaopt=area;
        r2opt=r2;
    end
    
    % ******************************************************************
    % test for convergence within maximum number of iterations allowed
    % ******************************************************************
    
    if abs(change) < evalin('base','mcr_als.alsOptions.opt.tolsigma'),
        
        %  finish the iterative optimization because convergence is achieved
        
%         disp(' ');disp(' ');
%         disp('CONVERGENCE IS ACHIEVED !!!!')
%         disp(' ')
%         disp(['Fitting error (lack of fit, lof) in % at the optimum = ', num2str(sdopt(1,1)),'(PCA) ', num2str(sdopt(1,2)), '(exp)']);
%         disp(['Percent of variance explained (r2)at the optimum is ',num2str(100*r2opt)]);
%         disp('Relative species conc. areas respect matrix (sample) 1at the optimum'),disp(rtopt')
%         disp(['Plots are at optimum in the iteration ', num2str(itopt)]);
        
        
        if isempty(evalin('base','mcr_als.alsOptions.out.out_conc'))==0
            assignin('base',evalin('base','mcr_als.alsOptions.out.out_conc'),copt)
        end
        if isempty(evalin('base','mcr_als.alsOptions.out.out_spec'))==0
            assignin('base',evalin('base','mcr_als.alsOptions.out.out_spec'),sopt)
        end
        if isempty(evalin('base','mcr_als.alsOptions.out.out_res'))==0
            assignin('base',evalin('base','mcr_als.alsOptions.out.out_res'),ropt)
        end
        if isempty(evalin('base','mcr_als.alsOptions.out.out_std'))==0
            assignin('base',evalin('base','mcr_als.alsOptions.out.out_std'),sdopt)
        end
        if isempty(evalin('base','mcr_als.alsOptions.out.out_area'))==0
            assignin('base',evalin('base','mcr_als.alsOptions.out.out_area'),areaopt)
        end
        if isempty(evalin('base','mcr_als.alsOptions.out.out_rat'))==0
            assignin('base',evalin('base','mcr_als.alsOptions.out.out_rat'),rtopt)
        end
        
        % save parameters
%         assignin('base','plot_optim_lof',sdopt(1,2));
%         assignin('base','plot_optim_R2',100*r2opt);
%         assignin('base','optim_niter',itopt);
%         assignin('base','optim_concs',copt);
%         assignin('base','optim_specs',sopt);
        
        
        evalin('base','mcr_als.alsOptions.resultats.plot_optim_lof=plot_optim_lof;');
        evalin('base','mcr_als.alsOptions.resultats.plot_optim_R2=plot_optim_R2;');
        evalin('base','mcr_als.alsOptions.resultats.optim_niter=optim_niter;');
        evalin('base','mcr_als.alsOptions.resultats.optim_concs=optim_concs;');
        evalin('base','mcr_als.alsOptions.resultats.optim_specs=optim_specs;');
        
        evalin('base','clear plot_optim_lof plot_optim_R2 optim_niter optim_concs optim_specs');
        
        
        als_end=1;
%         assignin('base','als_end',als_end);
%         assignin('base','copt_xxx',copt);
%         assignin('base','sopt_xxx',sopt);
%         assignin('base','sdopt_xxx',sdopt);
%         assignin('base','r2opt_xxx',r2opt);
%         assignin('base','rtopt_xxx',rtopt);
%         assignin('base','itopt_xxx',itopt);
%         assignin('base','change_xxx',sigman); % for std dev res vs exp
        als_res;
        evalin('base','clear als_end copt_xxx sopt_xxx sdopt_xxx r2opt_xxx rtopt_xxx itopt_xxx change_xxx');
        
        %         appCorrelation=evalin('base','mcr_als.alsOptions.correlation.appCorrelation;');
        %         if appCorrelation==1;
        %             assignin('base','correl_yout',yout);
        %             assignin('base','correl_ycal',ycal);
        %             assignin('base','correl_stats',stats);
        %             evalin('base','mcr_als.alsOptions.resultats.correlation.yout=correl_yout;');
        %             evalin('base','mcr_als.alsOptions.resultats.correlation.ycal=correl_ycal;');
        %             evalin('base','mcr_als.alsOptions.resultats.correlation.stats=correl_stats;');
        %             evalin('base','clear correl_yout correl_ycal correl_stats');
        %         end
        
        return         % 1st return (end of the optimization, convergence)
    end
    
    %  finish the iterative optimization if divergence occurs 20 times consecutively
    if idev > 50,
        %     if idev > 20,  original
%         disp(' ');disp(' ');
%         disp('FIT NOT IMPROVING FOR 50 TMES CONSECUTIVELY (DIVERGENCE?), STOP!!!')
%         disp(' ')
%         disp(['Fitting error (lack of fit, lof) in % at the optimum = ', num2str(sdopt(1,1)),'(PCA) ', num2str(sdopt(1,2)), '(exp)']);
%         disp(['Percent of variance explained (r2)at the optimum is ',num2str(100*r2opt)]);
%         disp('Relative species conc. areas respect matrix (sample) 1 at the optimum'),disp(rtopt)
%         disp(['Plots are at optimum in the iteration ', num2str(itopt)]);
        
        
%         if isempty(evalin('base','mcr_als.alsOptions.out.out_conc'))==0
%             assignin('base',evalin('base','mcr_als.alsOptions.out.out_conc'),copt)
%         end
%         if isempty(evalin('base','mcr_als.alsOptions.out.out_spec'))==0
%             assignin('base',evalin('base','mcr_als.alsOptions.out.out_spec'),sopt)
%         end
%         if isempty(evalin('base','mcr_als.alsOptions.out.out_res'))==0
%             assignin('base',evalin('base','mcr_als.alsOptions.out.out_res'),ropt)
%         end
%         if isempty(evalin('base','mcr_als.alsOptions.out.out_std'))==0
%             assignin('base',evalin('base','mcr_als.alsOptions.out.out_std'),sdopt)
%         end
%         if isempty(evalin('base','mcr_als.alsOptions.out.out_area'))==0
%             assignin('base',evalin('base','mcr_als.alsOptions.out.out_area'),areaopt)
%         end
%         if isempty(evalin('base','mcr_als.alsOptions.out.out_rat'))==0
%             assignin('base',evalin('base','mcr_als.alsOptions.out.out_rat'),rtopt)
%         end
        
        % save parameters
%         assignin('base','plot_optim_lof',sdopt(1,2));
%         assignin('base','plot_optim_R2',100*r2opt);
%         assignin('base','optim_niter',itopt);
%         assignin('base','optim_concs',copt);
%         assignin('base','optim_specs',sopt);
        
        evalin('base','mcr_als.alsOptions.resultats.plot_optim_lof=plot_optim_lof;');
        evalin('base','mcr_als.alsOptions.resultats.plot_optim_R2=plot_optim_R2;');
        evalin('base','mcr_als.alsOptions.resultats.optim_niter=optim_niter;');
        evalin('base','mcr_als.alsOptions.resultats.optim_concs=optim_concs;');
        evalin('base','mcr_als.alsOptions.resultats.optim_specs=optim_specs;');
        
        evalin('base','clear plot_optim_lof plot_optim_R2 optim_niter optim_concs optim_specs');
        
        als_end=2;
%         assignin('base','als_end',als_end);
%         assignin('base','copt_xxx',copt);
%         assignin('base','sopt_xxx',sopt);
%         assignin('base','sdopt_xxx',sdopt);
%         assignin('base','r2opt_xxx',r2opt);
%         assignin('base','rtopt_xxx',rtopt);
%         assignin('base','itopt_xxx',itopt);
%         assignin('base','change_xxx',sigman); % for std dev res vs exp
        als_res;
        evalin('base','clear als_end copt_xxx sopt_xxx sdopt_xxx r2opt_xxx rtopt_xxx itopt_xxx change_xxx');
        
        
        %         appCorrelation=evalin('base','mcr_als.alsOptions.correlation.appCorrelation;');
        %         if appCorrelation==1;
        %             assignin('base','correl_yout',yout);
        %             assignin('base','correl_ycal',ycal);
        %             assignin('base','correl_stats',stats);
        %             evalin('base','mcr_als.alsOptions.resultats.correlation.yout=correl_yout;');
        %             evalin('base','mcr_als.alsOptions.resultats.correlation.ycal=correl_ycal;');
        %             evalin('base','mcr_als.alsOptions.resultats.correlation.stats=correl_stats;');
        %             evalin('base','clear correl_yout correl_ycal correl_stats');
        %         end
        
        return          % 2nd return (end of optimization, divergence)
        
    end
    
    % this end refers to number of iterations initially proposed exceeded
    


% finish the iterative optimization if maximum number of allowed iterations is exceeded

%disp(' ');disp(' ');
% disp('NUMBER OF ITERATIONS EXCEEDED THE ALLOWED!')
% disp(' ')
% disp(['Fitting error (lack of fit, lof) in % at the optimum = ', num2str(sdopt(1,1)),'(PCA) ', num2str(sdopt(1,2)), '(exp)']);
% disp(['Percent of variance explained (r2)at the optimum is ',num2str(100*r2opt)]);
% disp('Relative species conc. areas respect matrix (sample) 1 at the optimum'),disp(rtopt)
% disp(['Plots are at optimum in the iteration ', num2str(itopt)]);
shiftSig = sopt;
shiftConc = copt;
% sub_conc = mean(copt,2);
% sub_abss = mean(sopt);
sub_conc = sum(copt,2);
sub_abss = sum(sopt);
% if isempty(evalin('base','mcr_als.alsOptions.out.out_conc'))==0
%     assignin('base',evalin('base','mcr_als.alsOptions.out.out_conc'),copt)
% end
% if isempty(evalin('base','mcr_als.alsOptions.out.out_spec'))==0
%     assignin('base',evalin('base','mcr_als.alsOptions.out.out_spec'),sopt)
% end
% if isempty(evalin('base','mcr_als.alsOptions.out.out_res'))==0
%     assignin('base',evalin('base','mcr_als.alsOptions.out.out_res'),ropt)
% end
% if isempty(evalin('base','mcr_als.alsOptions.out.out_std'))==0
%     assignin('base',evalin('base','mcr_als.alsOptions.out.out_std'),sdopt)
% end
% if isempty(evalin('base','mcr_als.alsOptions.out.out_area'))==0
%     assignin('base',evalin('base','mcr_als.alsOptions.out.out_area'),areaopt)
% end
% if isempty(evalin('base','mcr_als.alsOptions.out.out_rat'))==0
%     assignin('base',evalin('base','mcr_als.alsOptions.out.out_rat'),rtopt)
% end

% save parameters
% assignin('base','plot_optim_lof',sdopt(1,2));
% assignin('base','plot_optim_R2',100*r2opt);
% assignin('base','optim_niter',itopt);
% assignin('base','optim_concs',copt);
% assignin('base','optim_specs',sopt);
% 
% evalin('base','mcr_als.alsOptions.resultats.plot_optim_lof=plot_optim_lof;');
% evalin('base','mcr_als.alsOptions.resultats.plot_optim_R2=plot_optim_R2;');
% evalin('base','mcr_als.alsOptions.resultats.optim_niter=optim_niter;');
% evalin('base','mcr_als.alsOptions.resultats.optim_concs=optim_concs;');
% evalin('base','mcr_als.alsOptions.resultats.optim_specs=optim_specs;');
% 
% evalin('base','clear plot_optim_lof plot_optim_R2 optim_niter optim_concs optim_specs');
% 
% als_end=3;
% assignin('base','als_end',als_end);
% assignin('base','copt_xxx',copt);
% assignin('base','sopt_xxx',sopt);
% assignin('base','sdopt_xxx',sdopt);
% assignin('base','r2opt_xxx',r2opt);
% assignin('base','rtopt_xxx',rtopt);
% assignin('base','itopt_xxx',itopt);
% assignin('base','change_xxx',sigman); % for std dev res vs exp
% als_res;
% evalin('base','clear als_end copt_xxx sopt_xxx sdopt_xxx r2opt_xxx rtopt_xxx itopt_xxx change_xxx');

% appCorrelation=evalin('base','mcr_als.alsOptions.correlation.appCorrelation;');
% if appCorrelation==1;
%     assignin('base','correl_yout',yout);
%     assignin('base','correl_ycal',ycal);
%     assignin('base','correl_stats',stats);
%     evalin('base','mcr_als.alsOptions.resultats.correlation.yout=correl_yout;');
%     evalin('base','mcr_als.alsOptions.resultats.correlation.ycal=correl_ycal;');
%     evalin('base','mcr_als.alsOptions.resultats.correlation.stats=correl_stats;');
%     evalin('base','clear correl_yout correl_ycal correl_stats');
% end
% if exist('funcname')==0
%     funcname='alsOptimization_fcm2_re';
% else
%     funcname(length(funcname)+1)='alsOptimization_fcm2_re';
% end
return      %3rd return (end of optimization, number of iterations exceeded)
end
