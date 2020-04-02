%------------------------ Header Comment -------------------------------%
% Author:       Jian Hou, Weitian Chen
% email:        wtchen@cuhk.edu.hk
% Version:      1.0
% Date:         2020-04-02 15:11:50
% Description:  a simple demo script for MPF-SL 
% Functions:    
%               1. calculate Rmpfsl and Delta Rmt for three tissues
% History: 
%               1. First built on 2020-04-02.
% Others: 
%---------------------------------------------------------------------%
close all; clear all; clc;
% parameter libs for living tissues: liver, cartilage, muscle and white matter
%       Ref.: Stanisz GJ, Odrobina EE, Pun J, Escaravage M, Graham SJ, Bronskill MJ, et al.: 
%                   <<T1, T2 relaxation and magnetization transfer in tissue at 3T>>. Magnetic Resonance in Medicine. 2005;54(3):507–12.  
lib_tiss =  {'LIVER',   'CARTILAGE',    'MUSCLE',   'WHITE MATTER'};
%tiss pars: liver       cartilage       muscle      white matter     
lib_R1a =   [1./0.812,  1./1.168,       1./1.412,   1./1.084];
lib_R2a =   [1./42,     1./27,          1./50,      1./69] .* 1e3;
lib_R1c =   lib_R1a;
lib_R2c =   [1./7.7,    1./8.3,         1./8.7,     1./10] .* 1e6;
lib_fc  =   [6.9,       17.1,           7.4,        13.9] ./ 100;
lib_kca =   [51,        57,             66,         23];

% select tissues we want to investigate
var_tiss = {'LIVER', 'CARTILAGE', 'WHITE MATTER'};
nvar_tiss = length(var_tiss);

% sequence parameters
fsl1_a = [100]; fo1_a = [1000]; Nall = [4];

for iitiss = 1:nvar_tiss
    % load tissue pars
    tiss_save = char(var_tiss(iitiss));
    iitissue = find(1 == strcmp(lib_tiss, tiss_save));
    
    if ~isempty(iitissue)
        R1a = lib_R1a(iitissue);
        R2a = lib_R2a(iitissue);
        R1b = lib_R1c(iitissue);
        R2b = lib_R2c(iitissue);
        T2b = 1./R2b;
        fb  = lib_fc(iitissue);
        kba = lib_kca(iitissue);
    else
        error('no specific tissue, pls check!');
        return;
    end
    
    kab=kba*fb; 
    disp(tiss_save); 
    
        
    % load sequence pars
    for inn = 1:length(Nall)
        N = Nall(inn);
        disp(sprintf('      N is %d', N)); 
        for ii1 = 1:length(fsl1_a)
            for jj1 = 1:length(fo1_a)
                fsl1 = fsl1_a(ii1); 
                fo1  = fo1_a(jj1); 
                % FSL2=FSL*N, FO2=FO*N, as the papaer claimed.
                fsl2=fsl1*N; 
                fo2= fo1*N; 

                theta = atan2(fsl1, fo1);      
                % Reff, R1rho component due to free water pool
                Reff = R1a*(cos(theta))^2 + R2a*(sin(theta))^2; 

                % The Super-Lorentzian lineshape function, can be replaced by Gaussian lineshape if the objective are solids or gels.
                %       Ref.: Morrison C, Mark Henkelman R.: 
                %                   << A Model for Magnetization Transfer in Tissues>>. Magn Reson Med. 1995;33(4):475-482. doi:10.1002/mrm.1910330404
                fun = @(u, T2b, deltaw) sqrt(2/pi).* exp(-2* (deltaw*T2b./(3*u.^2-1)).^2)./(abs(3*u.^2-1)); 

                % Rrfc, the saturation rate term.
                Rrfc1 = integral(@(x)fun(x,T2b, 2*pi*fo1), 0, 1);
                Rrfc11 = (2*pi*fsl1)^2*pi*Rrfc1*T2b; 

                Rrfc2 = integral(@(x)fun(x,T2b, 2*pi*fo2), 0, 1);
                Rrfc22 = (2*pi*fsl2)^2*pi*Rrfc2*T2b; 

                % r1b1, r1b2, r1a and r2a, preparation for next-step calculation.
                r1b1 = R1b - Reff + Rrfc11;
                r1b2 = R1b - Reff + Rrfc22;

                r1a = R1a - Reff; 
                r2a = R2a - Reff; 

                for icase = 1:2
                    if icase==1
                        r1b = r1b1; 
                        deltaw = 2*pi*fo1;
                        w1 = 2*pi*fsl1; 
                    else
                        r1b = r1b2; 
                        deltaw = 2*pi*fo2;
                        w1 = 2*pi*fsl2; 
                    end

                    % Calculation for Rmt
                    %       Ref.: Zaiss M, Zu Z, Xu J, Schuenke P, Gochberg DF, Gore JC, et al.: 
                    %               <<A combined analytical solution for chemical exchange saturation transfer and semi-solid magnetization transfer>>. NMR in Biomedicine. 2015;28(2):217–30. 
                    % num, numerator of Rmt
                    num = (deltaw^2+r2a^2)*(kba*r1a + r1b*(kab+r1a)) + w1^2*r2a*(kba+r1b); 
                    % deno, denominator of Rmt
                    deno = (deltaw^2+r2a^2)*(kab+kba+r1a+r1b) + 2*r2a*(kba*r1a+r1b*(kab+r1a)) + w1^2*(r2a+kba+r1b); 
                    % Rmt for SL1 and SL2
                    rmt(icase) = num/deno; 
                end
                
                % d_rmt_full, The exact difference of Rmt
                d_rmt_full(ii1,jj1,inn) = rmt(2) - rmt(1); 

                % Our proposed Rmpfsl: d_rmt_app, The approximate difference of Rmt
                term1 = (1+fb)*kba + Rrfc11; 
                term1 = 1./term1; 
                term2 = (1+fb)*kba + Rrfc22; 
                term2 = 1./term2; 
                d_rmt_app(ii1,jj1,inn) = kba^2*fb*(1+fb)*(term1-term2); 
            end
        end
        d_rmt_full_inn = d_rmt_full(:,:,inn);
        d_rmt_app_inn  = d_rmt_app(:,:,inn);
    end
    save(['res_demo_MPFSL_', tiss_save, '.mat'], ...
            'fo1_a', 'fsl1_a', 'd_rmt_full', 'd_rmt_app',...
            'tiss_save', 'Nall');
end