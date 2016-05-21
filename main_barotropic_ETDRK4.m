clear all;
clc

global del2 ksq KX KY Nx Ny beta mu nu4 nu2 dpsidx k0 FILTER;


Nx=2*128;Ny=Nx;
Lx=2*pi;Ly=2*pi;
cl = get(gca,'colororder');

figNo=10; np=round(Nx/(64/3));

 mu = 0e-2;
nu4 = 0e-07; % e.g. N=128 choose 1e-06
nu2 = 0*.0005;
nofilter=0; % if nofilter=0 then it puts FILTER, if nofilter=1 it does not put FITLER
beta= 0;

bstr  =num2str(beta,'%1.2f');   bstr  =strrep(bstr,'.','p');
mustr =num2str(mu,'%1.2e');     mustr =strrep(mustr,'.','p');   mustr =strrep(mustr,'+','');    mustr =strrep(mustr,'-','m');
nu2str=num2str(nu2,'%1.2e');    nu2str=strrep(nu2str,'.','p');  nu2str=strrep(nu2str,'+','');   nu2str=strrep(nu2str,'-','m');
filenamesav = ['NL_barotr_b_' bstr '_mu_' mustr '_nu2_' nu2str '_N' num2str(Ny) '.mat'];




zmax = 40; % for plot
psimax = 1.5; % for plot


dt=0.005;
tfin=40;

Nthov=round(5/dt);
NtPsit=round(20/dt);
Nplot=round(1/dt);



dx=Lx/Nx;x=0:dx:Lx-dx;x=x.';
dy=Ly/Ny;y=0:dy:Ly-dy;y=y.';

kx=2*pi/Lx*[0:Nx/2-1 -Nx/2:-1];kx=kx.';
ky=2*pi/Ly*[0:Ny/2-1 -Ny/2:-1];ky=ky.';

[ X, Y]=meshgrid( x, y);
[KX,KY]=meshgrid(kx,ky);
ksq=KX.^2+KY.^2;
del2=-ksq;
del2(KX==0&KY==0)=1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTER

% a is chosen so that the energy at the largest nondim
% wavenumber K*dx be zero whithin machine double precision
s=4;
Kmax = Ny/2; Kmax_s=Kmax*dy;
kcut = 2/3*Kmax;kcut_s=kcut*dy;
a = -log(1e-15)/(Kmax_s-kcut_s)^s * dy^s;
K=sqrt(KX.^2+KY.^2);

FILTER = 1*ones(Ny,Nx).*abs(K<=Ny/3) + exp(-a*(K-kcut).^s).*abs(K>Ny/3);

if nofilter==1, FILTER = ones(Ny,Nx);end;
FILTERy = FILTER(:,1);

% figure(1); hold on;plot(fftshift(kx),fftshift(FILTER(1,:)),'--r','linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



T = 0:dt:tfin;
Ep=0*T;Ut=0*T;
Z =0*T;

UT   = zeros(Ny,round(length(T)/Nthov));
PsiT = zeros(Ny*Nx,round(length(T)/NtPsit+1));

% initial condition
Wkl=randn(size(KX))+1i*randn(size(KX));
Wkl = Wkl.*FILTER; Wkl(KX==0&KY==0)=0;

% zeta=exp(-((X-pi/2).^2+(Y-pi/2).^2)/.5^2) - exp(-((X-3*pi/2).^2+(Y-3*pi/2).^2)/.5^2);
% zeta = zeta-mean(zeta(:));
% zhat = fft2(zeta);

% this is mcwilliams initial condition |\hat\psi|^2
psikl = sqrt( 1./sqrt(ksq) .* 1./(1+(ksq.^2/6^4)) ).*FILTER;
psikl(1,1)=0;

zhat = del2.*(randn(Ny,Nx)+1i*randn(Ny,Nx)).*psikl;
% zeta = zeta - mean(zeta(:));
zeta = real(ifft2(zhat));
zhat = fft2(zeta);


psihat = zhat./del2;
u = real(ifft2(-1i*KY.*psihat));
v = real(ifft2(+1i*KX.*psihat));


E0=.5*sum(u(:).^2+v(:).^2)*dx*dy/(Lx*Ly);

zhat = zhat*sqrt(.5/E0);
psihat = zhat./del2;
u = real(ifft2(-1i*KY.*psihat));
v = real(ifft2(+1i*KX.*psihat));


Ep(1)= -.5*sum(abs(zhat(:)).^2./del2(:))/(Nx*Ny)^2;
Z(1) =  .5*sum(abs(zhat(:)).^2)/(Nx*Ny)^2;

%%

ithov = 1;itpsi=1;
Thov(ithov) = T(1);
UT(:,ithov) = mean(u,2);
ithov = ithov+1;


teddy = 1/sqrt(Z(1));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of the EDTRK4 coefficients

L = -1i*beta*KX./del2 - mu - nu2*ksq - nu4*ksq.^2;
eL  = exp(L*dt);
eL2 = exp(L*dt/2);
% IL = 1./L;
% ILeL2 = IL.*(eL2-1);


a1=1;
a2=1;
a3=1;

M=2*64;
% r = ceil(max(abs(L(:))))*exp(1i*pi*((1:M)-.5)/M);
% r = 1*exp(1i*pi*((1:M)-.5)/M);
r = 1*exp(2i*pi*((1:M)/M));
E = exp(L*dt); E2 = exp(L*dt/2);
fu = zeros(Ny,Nx); fab = fu; fc = fu; Q = fu;
tic
for j = 1:M
     z = r(j) + L*dt;
     Q =   Q + dt*( exp(z/2)-1 )./z;
    fu =  fu + dt*( -4 -z +exp(z).*(4 -3*z +z.^2) )./z.^3;
   fab = fab + dt*( +2 +z +exp(z).*(-2 +z) )./z.^3;
    fc =  fc + dt*( -4 -3*z -z.^2 +exp(z).*(4 -z) )./z.^3;
end
% f1 = real(f1/M); f2 = real(f2/M); f3 = real(f3/M); Q = real(Q/M);
fu = (fu/M); fab = (fab/M); fc = (fc/M); Q = (Q/M);
% format long
% [fu(4,5);fab(4,5);fc(4,5);Q(4,5)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


psi = real(ifft2(zhat./del2));
PsiT(:,itpsi) = reshape(psi,Nx*Ny,1);
itpsi=itpsi+1;

%% plot initial conditions

it=1;

psihat = zhat./del2;
    zeta=real(ifft2(zhat));
    psi = real(ifft2(psihat));
    u = real(ifft2(-1i*KY.*psihat));
    v = real(ifft2(+1i*KX.*psihat));
    cfl = sqrt(max((u(:)).^2+v(:).^2))*dt/dx;
    display(['t=' num2str(T(it)) '  cfl=' num2str(cfl,'%1.3f') '  E=' num2str(Ep(1),'%1.9f') ]);

    nRow=2;nCol=2;

    figure(figNo);clf;
    subplot(nRow,nCol,2)
    pcolor(X,Y,zeta);shading interp;
    axis square;colorbar;
	%caxis([-1 1]*zmax)
    axis([0 x(end) 0 y(end)]);
    title(['$\zeta(\mathbf{x},t=' num2str(T(it),'%1.1f') ')$'],'fontsize',18,'interpreter','latex');        xlabel('$x$','fontsize',18,'interpreter','latex');
    ylabel('$y$','fontsize',18,'interpreter','latex');
    set(gca,'tickLength',[.025 .025],'xtick',[0:2:10],'ytick',[0:2:10],'layer','top')

    subplot(nRow,nCol,1)
    pcolor(X,Y,psi);shading interp;
    hold on;
    quiver(X(1:np:end,1:np:end),Y(1:np:end,1:np:end),u(1:np:end,1:np:end),v(1:np:end,1:np:end),'k','linewidth',1)
    hold off;
    %caxis([-1 1]*psimax)
    axis square;colorbar;
    axis([0 x(end) 0 y(end)]);
    title(['$\psi(\mathbf{x},t=' num2str(T(it),'%1.1f') ')$'],'fontsize',18,'interpreter','latex');
    xlabel('$x$','fontsize',18,'interpreter','latex');
    ylabel('$y$','fontsize',18,'interpreter','latex');
    set(gca,'tickLength',[.025 .025],'xtick',[0:2:10],'ytick',[0:2:10],'layer','top')
    rgb1=[0 0 .8];
    rgb2=[.8 0 0];
    s=linspace(0,1,151);
    cmap = diverging_map(s,rgb1,rgb2);
    colormap(cmap);

%         subplot(nRow,nCol,[5,6])
% %         contourf(Thov(1:ithov-1),y,UT(:,1:ithov-1),'linestyle','none');colorbar;
%         pcolor(mu*Thov(1:1:ithov-1),y,UT(:,1:1:ithov-1));shading interp;colorbar;
%         %caxis([-1 1]*max(max(abs(UT(:,1:ithov-1)))));
%         xlabel('$\mu t$','fontsize',18,'interpreter','latex');
%         ylabel('$U+\overline{u}$','fontsize',18,'interpreter','latex');
%         xlim(mu*[0 max(Thov(1:1:ithov-1))]);



    drawnow;

%%


EE = -.5*abs(zhat).^2./del2/(Nx*Ny)^2;
EE = +.5*abs(zhat).^2/(Nx*Ny)^2;


dkx = kx(2)-kx(1);
dky = ky(2)-ky(1);
dkr = sqrt(dkx^2+dky^2);
kmax=min(max(kx(:)),max(ky(:)));
Kr = dkr/2:dkr:kmax;
wv = sqrt(KX.^2+KY.^2);



tic
for it=2:length(T);


    [nlin0_z] = Dtzeta_EDTRK4(zhat);
    k1z = eL2.*zhat + Q.*nlin0_z;

    [nlin1_z] = Dtzeta_EDTRK4(k1z);
    k2z = eL2.*zhat + Q.*nlin1_z;

    [nlin2_z] = Dtzeta_EDTRK4(k2z);
    k3z = eL2.*k1z + Q.*(2*nlin2_z-nlin0_z);

    [nlin3_z] = Dtzeta_EDTRK4(k3z);


    zhatnew  = eL.*zhat + fu.*nlin0_z + 2*fab.*(nlin1_z+nlin2_z) +fc.*nlin3_z ;
    zhat = zhatnew.*FILTER;

    Ep(it) = -.5*sum(abs(zhat(:)).^2./del2(:))/(Nx*Ny)^2;
    Z(it)  =  .5*sum(abs(zhat(:)).^2)/(Nx*Ny)^2;

    EE = EE*(1-1/it) -.5*abs(zhat).^2./del2/(Nx*Ny)^2/it;
    EE = EE*(1-1/it) +.5*abs(zhat).^2/(Nx*Ny)^2/it;

    if rem(it,Nthov)==1
        psihat = zhat./del2;
        zeta=real(ifft2(zhat));
        psi = real(ifft2(psihat));
        u = real(ifft2(-1i*KY.*psihat));
        v = real(ifft2(+1i*KX.*psihat));
        Thov(ithov) = T(it);
        UT(:,ithov) =  mean(u,2);
        ithov=ithov+1;
    end

    if rem(it,NtPsit)==1||it==length(T)
        psi = real(ifft2(zhat./del2));
        PsiT(:,itpsi) = reshape(psi,Nx*Ny,1);
        itpsi=itpsi+1;
    end

%     if rem(it,100*Nthov)==1
%         cfl = sqrt(max((U+u(:)).^2+v(:).^2))*dt/dx;
%         if cfl>.7
%             errormsg='cfl was larger than .7';
%             save(filenamesav);
%             error('cfl was larger than .7')
%         end
%         display(['t=' num2str(T(it)) '  cfl=' num2str(cfl,'%1.3f')]);
%     end

    if rem(it,Nplot)==1||it==length(T)

            Ekl=EE;
            Er=0*Kr;
            for ir = 1:length(Kr)
                kr = Kr(ir);
                fkr = ones(size(wv)).*(wv>=kr-dkr/2 & wv<= kr+dkr/2);
                fkr=round(fkr);
                dth = 2*pi/(sum(fkr(:))-1);
                Er_int = Ekl.*wv.^1.*fkr;
                Er(ir) = sum(Er_int(:))*dth;
            end
%%
            figure(298)
            loglog(Kr,Er,'*-','linewidth',3)

            hold on;
            plot(Kr,40*Kr.^(-5/3),'-k');
            hold off;
            xlim([1e0 Nx])

            ylim([1e-7 1e-1])
            ylim([1e-3 1e1])

        %%
        psihat = zhat./del2;
        zeta=real(ifft2(zhat));
        psi = real(ifft2(psihat));
        u = real(ifft2(-1i*KY.*psihat));
        v = real(ifft2(+1i*KX.*psihat));
        cfl = sqrt(max((u(:)).^2+v(:).^2))*dt/dx;
        display(['t=' num2str(T(it)) '  cfl=' num2str(cfl,'%1.3f') '  E=' num2str(Ep(it),'%1.9f') ]);

        nRow=2;nCol=2;

        figure(figNo);clf;
        subplot(nRow,nCol,2)
        pcolor(X,Y,zeta);shading interp;
        axis square;colorbar;
        caxis([-1 1]*zmax)
        axis([0 x(end) 0 y(end)]);
        title(['$\zeta(\mathbf{x},t=' num2str(T(it),'%1.1f') ')$'],'fontsize',18,'interpreter','latex');        xlabel('$x$','fontsize',18,'interpreter','latex');
        ylabel('$y$','fontsize',18,'interpreter','latex');
        set(gca,'tickLength',[.025 .025],'xtick',[0:2:10],'ytick',[0:2:10],'layer','top')
%         clim = get(gca,'CLim');cmax = min(abs(clim));caxis([-1 1]*cmax);
        subplot(nRow,nCol,1)
        pcolor(X,Y,psi);shading interp;
        hold on;
        quiver(X(1:np:end,1:np:end),Y(1:np:end,1:np:end),u(1:np:end,1:np:end),v(1:np:end,1:np:end),'k','linewidth',1)
        hold off;
%         caxis([-1 1]*psimax)
        axis square;colorbar;
        axis([0 x(end) 0 y(end)]);
        title(['$\psi(\mathbf{x},t=' num2str(T(it),'%1.1f') ')$'],'fontsize',18,'interpreter','latex');
        xlabel('$x$','fontsize',18,'interpreter','latex');
        ylabel('$y$','fontsize',18,'interpreter','latex');
        set(gca,'tickLength',[.025 .025],'xtick',[0:2:10],'ytick',[0:2:10],'layer','top')
        clim = get(gca,'CLim');cmax = min(abs(clim));caxis([-1 1]*cmax);
        subplot(nRow,nCol,[3,4])
        plot(T(1:it),Z(1:it)/Z(1),'-b','linewidth',2);
        hold on;
        plot(T(1:it),Ep(1:it)/Ep(1),'-r','linewidth',2);
        hold off;
        xlim([0 T(it)])
        title('b: $Z/Z(t=0)$~,~r: $E/E(t=0)$','fontsize',18,'interpreter','latex');
        xlabel('$t$','fontsize',18,'interpreter','latex');
        xlim([0 T(it)]);

        drawnow;



    end

%     if rem(it,round(length(T)/20))==1||it==length(T)
%         save(filenamesav);
%     end
end
toc

% save(filenamesav);
