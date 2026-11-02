!=============================================================================
!                   SUPERMEGAGIGAMODULECOOLINGQUIDEPOTE
!=============================================================================
! Les subroutines et fonctions d'interet general sont :
!
! ROUTINES A APPELER PAR LE CODE HYDRO
!
!    subroutine set_model(...) :
!          pour choisir le modele de cooling et ses parametres
!
!    subroutine set_table(aexp) :
!          pour creer la table avec les parametres par defaut
!          Plus pratique a appeler que
!          cmp_table(nH_min,nH_max,T2_min,T2_max,nbin_n,nbin_T,aexp)
!
!    subroutine solve_cooling(...) :
!          pour calculer le cooling
!
! ROUTINE A MODIFIER SI NECESSAIRE
!
!    function J0simple(aexp) :
!          donne le J0 en fonction du redshift dans les modeles Teyssier
!          ou Theuns
!
!=============================================================================
module cooling_module
  use constants
  implicit none
  logical,parameter :: verbose_cooling=.false.
  real(kind=8),parameter ::smallnum_cooling= 1d-30
  integer,parameter :: HI      = 1
  integer,parameter :: HEI     = 2
  integer,parameter :: HEII    = 3
  ! Utilisation de table%n_spec si necessaire
  logical, parameter :: if_species_abundances=.true.
  ! Facteur correctif de Theuns et al.
  real(kind=8),parameter :: dumfac_ion_theuns=2d0
  real(kind=8),parameter :: dumfac_rec_theuns=0.75D0  ! idem
  ! Les parametres de la table par defaut
  integer,parameter     :: nbin_T_fix=101
  integer,parameter     :: nbin_n_fix=161
  real(kind=8),parameter:: nH_min_fix=1d-10
  real(kind=8),parameter:: nH_max_fix=1d+6
  real(kind=8),parameter:: T2_min_fix=1d-2
  real(kind=8),parameter:: T2_max_fix=1d+9
  ! Order polynomial interpolation
  integer,parameter :: Nordercourty=7

  type cooling_table
     integer::n1
     integer::n2
     real(kind=8),dimension(:)    ,allocatable::nH
     real(kind=8),dimension(:)    ,allocatable::T2
     real(kind=8),dimension(:,:)  ,allocatable::cool
     real(kind=8),dimension(:,:)  ,allocatable::heat
     real(kind=8),dimension(:,:)  ,allocatable::cool_com
     real(kind=8),dimension(:,:)  ,allocatable::heat_com
     real(kind=8),dimension(:,:)  ,allocatable::metal
     real(kind=8),dimension(:,:)  ,allocatable::cool_prime
     real(kind=8),dimension(:,:)  ,allocatable::heat_prime
     real(kind=8),dimension(:,:)  ,allocatable::cool_com_prime
     real(kind=8),dimension(:,:)  ,allocatable::heat_com_prime
     real(kind=8),dimension(:,:)  ,allocatable::metal_prime
     real(kind=8),dimension(:,:)  ,allocatable::mu
     real(kind=8),dimension(:,:,:),allocatable::n_spec
  end type cooling_table

  type cooling_t

     type(cooling_table)::table

     real(kind=8) :: X = 0.76d0
     real(kind=8) :: Y = 0.24d0
     real(kind=8) :: mu_mol = 1.2195d0
     real(kind=8) :: dumfac_ion=dumfac_ion_theuns
     real(kind=8) :: dumfac_rec=dumfac_rec_theuns

     logical :: first = .true.

     ! On DOIT AVOIR OU teyssier OU theuns OU madau
     ! OU courty avec un OU exclusif
     logical :: teyssier=.false.
     logical :: theuns=.false.
     logical :: madau=.false.
     logical :: courty=.true.  ! Default model

     ! Si teyssier ou theuns :
     real(kind=8) :: J0in=1d-22  ! J0 default
     real(kind=8) :: J0min=1d-29 ! Valeur minimale du J0
     logical :: force_j0_one=.false. ! Force constant UV bkg
     ! (saturation a grand redshift)
     real(kind=8) :: aexp_ref=0.0001d0
     real(kind=8) :: J0min_ref=2.77168510365299962D-25 ! J0min_ref precalcule pour
     ! H0=70, omegab=0.04, omega0=0.3, omegaL=0.7
     logical :: high_z_realistic_ne=.true. ! Calcul du J0min de telle sorte
     ! que le n_e soit realiste a grand z. J0min=J0min_ref/(aexp/aexp_ref)^2
     real(kind=8) :: alpha=1d0   ! J(nu) \propto nu^{-alpha}
     ! Si madau :
     real(kind=8) :: normfacJ0=0.74627d0   ! Facteur de normalisation pour J0
     ! pour un J(nu,z) de type haardt et Madau
     ! Ce facteur la est celui utilise par Dave et al. pour LCDM
     ! Sauvegarde des termes de cooling/heating dans les

     ! variables en dessous
     real(kind=8)::cb1s,cb2s,cb3s,ci1s,ci2s,ci3s,cr1s,cr2s,cr3s,cds
     real(kind=8)::ce1s,ce3s,ch1s,ch2s,ch3s,cocs,cohs
     real(kind=8)::cool_out, heat_out

     real(kind=8) :: zreioniz=8.5d0
     ! Ordre+1 du polynome d'interpolation (NE PAS CHANGER)
     real(kind=8) :: coefcourty(0:Nordercourty,6)= reshape(           &
          (/-13.5857d0,     1.24475d0,    0.187739d0,  &
          -0.430409d0,    0.152544d0,  -0.0246448d0, &
          0.00192622d0, -5.89772d-05,               &
          -14.0242d0,     1.99211d0,   -0.490766d0,  &
          -0.122646d0,    0.0776501d0, -0.0146310d0, &
          0.00123335d0, -3.96066d-05,               &
          -15.6627d0,     0.128240d0,   1.65633d0,   &
          -1.23799d0,     0.372157d0,  -0.0561687d0, &
          0.00422696d0, -0.000126344d0,             &
          -24.8422d0,     1.50750d0,   -0.0699428d0, &
          -0.308682d0,    0.122196d0,  -0.0205179d0, &
          0.00163695d0, -5.08050d-05,               &
          -25.0252d0,     1.79577d0,   -0.159054d0,  &
          -0.300924d0,    0.125343d0,  -0.0214598d0, &
          0.00173377d0, -5.43576d-05,               &
          -26.4168d0,     0.0479454d0,  1.70948d0,   &
          -1.26395d0,     0.378922d0,  -0.0570957d0, &
          0.00428897d0, -0.000127909d0 /),(/Nordercourty+1,6/) )

     ! UVB change tracking
     real(kind=8) :: uvb_delta = 0.05d0      ! Relative change threshold
     real(kind=8) :: uvb_last_update = -1.0d0 ! UVB HI rate when table was last computed

  end type cooling_t

contains
!=======================================================================
subroutine set_model(c,Nmodel,J0in_in,J0min_in,alpha_in,normfacJ0_in,  &
 &                   zreioniz_in, correct_cooling,realistic_ne,h,      &
 &                   omegab,omega0,omegaL,astart_sim,T2_sim,mu_mol)
!=======================================================================
! Nmodel(integer) =1 : Teyssier : ancien choix de l'evolution et de la forme du J(nu,z)
!                  2 : Theuns   : pareil mais avec les fonctions interpolees de Theuns (+ rapide)
!                  3 : Madau    : J(nu,z) de Theuns et al. 1998 avec les anciennes mesures de
!                                 Haardt et Madau (HM)
!                  6 : Courty
!                 -1 : defaut defini dans le module
! J0in_in (dble) : valeur du J0 utilisee pour Teyssier et Theuns
!            Exemple : J0in_in=1d-22
!            J0in_in <= 0 utilise le defaut defini dans le module
! J0min_in (dble) : valeur du J0min ou J0min_ref (voir option realistic_ne)
!            utilisee dans tous les modeles a grand redshift
!            Exemple : J0min_in=1d-29
!            J0min_in <= 0 utilise le defaut defini dans le module
! alpha_in (dble) : valeur de l'indice spectral du J(nu) \propto nu^{-alpha}
!            Exemple : alpha=1.
!            alpha_in < 0 utilise le defaut defini dans le module
! zreioniz_in (dble) : valeur du redshift de reionisation
!            Exemple : zerion=10.
!            zreioniz_in < 0 utilise le defaut defini dans le module
! normfacJ0_in (dble) : valeur du facteur de normalisation dans le cas des
!            spectres de Haardt et Madau. C'est un nombre de l'ordre de
!            l'unite en general plus petit que 1.
!            Exemple : normfacJ0_in=0.74627
!            normfacJ0_in prend le defaut defini dans le module
! correct_cooling (integer) : 0 : pas de correction
!                             1 : correction de Theuns et al 98
!                            -1 : defaut defini dans le module
! realistic_ne (integer) : 0 : pas de n_e realiste a grand redshift :
!                              Le J0min reste le meme quel que soit le redshift
!                              (J0min=J0min_in si celui-ci est > 0)
!                          1 : n_e realiste a grand redshift : J0min proportionnel a 1/a^2
!                              egal initialement a J0min_ref pour a=aexp_ref=0.0001
!                              (J0min_ref=J0min_in si celui-ci est > 0)
!                          2 : RECOMMANDE : pareil que 1, mais J0min_ref est calcule de
!                              maniere iterative pour avoir le bon n_e a z=19.
!                              Le J0min_in n'est pas relevant dans ce cas la.
! h (dble)          : H0/100
! omegab (dble)     : omega baryons
! omega0 (dble)     : omega matiere total
! omegaL (dble)     : omega Lambda
! astart_sim (dble) : redshift auquel on veut commencer la simulation
! T2_sim     (dble) : ce sera en output, le T/mu en K a ce redshift pour des regions de contraste
!                     de densite nul.
! mu_mol (dble)     : mean molecular weight
!
! NOTE :
! Dans les cas madau, le J0 a grand redshift est calcule comme
! dans l'option theuns :
!   madau :       pour z >= 15 ou quand le taux trouve est plus petit que celui donne par
!                 l'option theuns=.true.
!   courty :      voire routine correspondante
!=======================================================================
  implicit none
  type(cooling_t)::c
  real(kind=8) :: J0in_in,zreioniz_in,J0min_in,alpha_in,normfacJ0_in,astart_sim,T2_sim
  real(kind=8) :: J0min_ref_calc,h,omegab,omega0,omegaL,mu_mol
  integer :: Nmodel,correct_cooling,realistic_ne
  real(kind=8) :: astart,aend,dasura,T2end,mu,ne,minus1
  if (Nmodel /= -1) then
     c%teyssier=.false.
     c%theuns=.false.
     c%madau=.false.
     c%courty=.false.
     if (Nmodel==1) then
        c%teyssier=.true.
     elseif (Nmodel==2) then
        c%theuns=.true.
     elseif (Nmodel==3) then
        c%madau=.true.
     elseif (Nmodel==6) then
        c%courty=.true.
     else
        write(*,*) 'ERROR in set_model : wrong value of Nmodel'
        write(*,*) 'Nmodel =',Nmodel
        STOP
     endif
  endif
  if (J0in_in >= 0.0) c%J0in=J0in_in
  if (zreioniz_in >= 0.0) c%zreioniz=zreioniz_in
  if (alpha_in > 0.0) c%alpha=alpha_in
  if (normfacJ0_in > 0.0) c%normfacJ0=normfacJ0_in
  if (correct_cooling == 0) then
     c%dumfac_ion=1d0
     c%dumfac_rec=1d0
  elseif (correct_cooling == 1) then
     c%dumfac_ion=dumfac_ion_theuns
     c%dumfac_rec=dumfac_rec_theuns
  elseif (correct_cooling /= -1) then
     write(*,*) 'ERROR in set_model : wrong value of correct_cooling'
     write(*,*) 'correct_cooling =',correct_cooling
     STOP
  endif
  if (realistic_ne == 0) then
     astart=5d-4
     c%high_z_realistic_ne=.false.
     if (J0min_in > 0d0) c%J0min=J0min_in
  elseif (realistic_ne == 1) then
     astart=c%aexp_ref
     c%high_z_realistic_ne=.true.
     if (J0min_in > 0d0) c%J0min_ref=J0min_in
  elseif (realistic_ne == 2) then
     astart=c%aexp_ref
     c%high_z_realistic_ne=.true.
     call compute_J0min(c,h,omegab,omega0,omegaL,J0min_ref_calc)
     c%J0min_ref=J0min_ref_calc
  else
     write(*,*) 'ERROR in set_model : wrong value of realistic_ne'
     write(*,*) 'realistic_ne =',realistic_ne
     STOP
  endif
  if (astart_sim < astart) then
     write(*,*) 'ERROR in set_model : astart_sim is too small.'
     write(*,*) 'astart     =',astart
     write(*,*) 'astart_sim =',astart_sim
     STOP
  endif
  c%mu_mol = mu_mol
  ! Calcul de la temperature initiale
  aend=astart_sim
  dasura=0.02d0
  minus1=-1
  call evol_single_cell(c,astart,aend,dasura,h,omegab,omega0,omegaL,minus1,T2end,mu,ne,.false.)
  if (verbose_cooling) write(*,*) 'Starting temperature in K :',T2end*mu
  T2_sim=T2end
end subroutine set_model
!=======================================================================
subroutine set_table(c,aexp)
!=======================================================================
  implicit none
  type(cooling_t)::c
  real(kind=8) :: aexp
  integer :: nbin_n,nbin_T
  real(kind=8) :: nH_min,nH_max,T2_min,T2_max
  nH_min=nH_min_fix
  nH_max=nH_max_fix
  T2_min=T2_min_fix
  T2_max=T2_max_fix
  nbin_n=nbin_n_fix
  nbin_T=nbin_T_fix
  call cmp_table(c,nH_min,nH_max,T2_min,T2_max,nbin_n,nbin_T,aexp)
end subroutine set_table
!=======================================================================
subroutine output_cool(c,filename)
!=======================================================================
  implicit none
  type(cooling_t)::c
  character(LEN=80)::filename
  open(unit=10,file=filename,form='unformatted')
  write(10)c%table%n1,c%table%n2
  write(10)c%table%nH
  write(10)c%table%T2
  write(10)c%table%cool
  write(10)c%table%heat
  write(10)c%table%cool_com
  write(10)c%table%heat_com
  write(10)c%table%metal
  write(10)c%table%cool_prime
  write(10)c%table%heat_prime
  write(10)c%table%cool_com_prime
  write(10)c%table%heat_com_prime
  write(10)c%table%metal_prime
  write(10)c%table%mu
  if (if_species_abundances) write(10)c%table%n_spec
  close(10)
end subroutine output_cool
!=======================================================================
subroutine evol_single_cell(c,astart,aend,dasura,h,omegab,omega0,omegaL, &
 &                          J0min_in,T2end,mu,ne,if_write_result)
!=======================================================================
! astart : valeur du facteur d'expansion au debut du calcul
! aend   : valeur du facteur d'expansion a la fin du calcul
! dasura : la valeur de da/a entre 2 pas de temps
! h      : la valeur de H0/100
! omegab : la valeur de Omega baryons
! omega0 : la valeur de Omega matiere (total)
! omegaL : la valeur de Omega Lambda
! J0min_in : la valeur du J0min a injecter :
!          Si high_z_realistic_ne alors c'est J0min a a=astart qui
!          est considere
!          Sinon, c'est le J0min habituel.
!          Si J0min_in <=0, les parametres par defaut ou predefinis
!          auparavant sont pris pour le J0min.
! T2end  : Le T/mu en output
! mu     : le poids moleculaire en output
! ne     : le ne en output
! if_write_result : .true. pour ecrire l'evolution de la temperature
!          et de n_e sur l'ecran.
!=======================================================================
  implicit none
  type(cooling_t)::c
  real(kind=8)::astart,aend,T2end,h,omegab,omega0,omegaL,J0min_in,ne,dasura
  logical :: if_write_result
  real(kind=8)::aexp,daexp,dt_cool,coeff,coeff2
  real(kind=8)::T2_com,T2_old,T2,T2_left,T2_right,err_T2
  real(kind=8)::nH_com,nH
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec
  real(kind=8) ::mu
  real(kind=8) ::cool_tot,heat_tot,cool_com,heat_com
  real(kind=8) ::diff
  integer::niter
  real(kind=8) :: n_spec(1:6)
  if (J0min_in > 0.0) then
     if (c%high_z_realistic_ne) then
        c%J0min_ref = J0min_in
        c%aexp_ref = astart
     else
        c%J0min = J0min_in
     endif
  endif
  aexp = astart
  T2_com = 2.726d0 / aexp * aexp**2 / c%mu_mol
  nH_com = omegab*rhoc*h**2*c%X/mH
  do while (aexp < aend)
     daexp = dasura*aexp
     dt_cool=daexp/(aexp*100*h*3.2408608d-20*HsurH0(1/aexp-1d0,omega0,omegaL,1d0-omega0-omegaL))

     nH = nH_com/aexp**3
     T2_old = T2_com/aexp**2

     ! Compute radiative ionization and heating rates
     call set_rates(c,t_rad_spec,h_rad_spec,aexp)

     ! Iteration to find new T2
     err_T2=1d0
     T2_left=1d-2
     T2_right=1d8
     niter=0
     coeff = 2d0*nH*c%X/3d0/kB
     coeff2 = 2d0*c%X/3d0/kB
     do while (err_T2 > 1d-10.and.niter <= 100)
        T2=0.5d0*(T2_left+T2_right)
        call cmp_cooling(c,T2,nH,t_rad_spec,h_rad_spec,cool_tot,heat_tot,cool_com,heat_com,mu,aexp,n_spec)
        diff = coeff*(heat_tot-cool_tot) + coeff2*(heat_com-cool_com) + (T2_old-T2)/dt_cool
        if(diff>0.)then
           T2_left =0.5d0*(T2_left+T2_right)
           T2_right=T2_right
        else
           T2_left =T2_left
           T2_right=0.5d0*(T2_left+T2_right)
        end if
        err_T2=abs(T2_right-T2_left)/T2_left
        niter=niter+1
     end do
     if (niter > 100) then
        write(*,*) 'ERROR in evol_single_cell : too many iterations'
        STOP
     endif
     T2_com=T2*aexp**2
     aexp = aexp + daexp
     if (if_write_result) write(*,'(4(1pe10.3))')aexp,nH,T2_com*mu/aexp**2,n_spec(1)/nH
  end do
  T2end=T2
  ne=n_spec(1)/nH
end subroutine evol_single_cell
!=======================================================================
subroutine compute_J0min(c,h,omegab,omega0,omegaL,J0min_in)
!=======================================================================
  implicit none
  type(cooling_t)::c
  real(kind=8) :: omega0,omegaL,h,omegab,ne_to_find,mu
  real(kind=8) :: astart,aend,J0min_in,T2end,ne
  real(kind=8) :: J0min_left,J0min_right,err_J0min,diff,xval,dasura
  integer :: niter
  logical :: if_write_result=.false.

  xval=sqrt(omega0)/(h*omegab)
  ne_to_find=1.2d-5*xval ! From the book of Peebles p. 173
  astart=c%aexp_ref
  aend=MIN(0.05d0,0.5d0/(1d0+c%zreioniz)) ! Always end before reionization
  dasura=0.05d0
  err_J0min=1d0
  J0min_left=1d-20
  J0min_right=1d-30
  niter=0
  do while (err_J0min > 1d-3 .and. niter <= 100)
     J0min_in=0.5d0*(J0min_left+J0min_right)
     call evol_single_cell(c,astart,aend,dasura,h,omegab,omega0,omegaL,J0min_in,T2end,mu,ne,if_write_result)
     diff=ne-ne_to_find
     if (diff>0d0) then
        J0min_left=0.5d0*(J0min_left+J0min_right)
        J0min_right=J0min_right
     else
        J0min_left=J0min_left
        J0min_right=0.5d0*(J0min_left+J0min_right)
     endif
     err_J0min=abs(J0min_right-J0min_left)/J0min_left
     niter=niter+1
  enddo
  if (niter > 100) then
     write(*,*) 'ERROR in compute_J0min : too many iterations'
     STOP
  endif
  if (verbose_cooling)  write(*,*) 'J0min found ',J0min_in
end subroutine compute_J0min
!=======================================================================
subroutine solve_cooling(c,nH,T2,zsolar,boost,dt,deltaT2,ncell)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer::ncell
  real(kind=8)::dt
  real(kind=8),dimension(1:ncell)::nH,T2,deltaT2,zsolar,boost

  real(kind=8)::facT,dlog_nH,dlog_T2,precoeff,h,h2,h3
  real(kind=8)::metal,cool,heat,cool_com,heat_com,yy,yy2,yy3
  real(kind=8)::metal_prime,cool_prime,heat_prime,cool_com_prime,heat_com_prime,wcool
  real(kind=8)::lambda,lambda_prime,logT2max
  real(kind=8)::fa,fb,fprimea,fprimeb,alpha,beta,gamma
  real(kind=8),dimension(1:ncell)::tau,tau_old
  real(kind=8),dimension(1:ncell)::time,time_old,facH,zzz,tau_ini
  real(kind=8),dimension(1:ncell)::w1H,w2H,wmax,time_max
  real(kind=8)::varmax=4d0
  integer::i,i_T2,iter,n,n_active
  integer,dimension(1:ncell)::ind,i_nH
  logical::tau_negative

  ! Initializations
  logT2max=log10(T2_max_fix)
  dlog_nH=dble(c%table%n1-1)/(c%table%nH(c%table%n1)-c%table%nH(1))
  dlog_T2=dble(c%table%n2-1)/(c%table%T2(c%table%n2)-c%table%T2(1))
  h=1d0/dlog_T2
  h2=h*h
  h3=h2*h
  precoeff=2d0*c%X/(3d0*kB)
  do i=1,ncell
     zzz(i)=zsolar(i)
     facH(i)=MIN(MAX(log10(nH(i)/boost(i)),c%table%nH(1)),c%table%nH(c%table%n1))
     i_nH(i)=MIN(MAX(int((facH(i)-c%table%nH(1))*dlog_nH)+1,1),c%table%n1-1)
     w1H(i)=(c%table%nH(i_nH(i)+1)-facH(i))*dlog_nH
     w2H(i)=(facH(i)-c%table%nH(i_nH(i)  ))*dlog_nH
     tau(i)=T2(i)
     tau_ini(i)=T2(i)
     time_max(i)=dt*precoeff*nH(i)
     time(i)=0d0
     wmax(i)=1d0/time_max(i)
     ind(i)=i
  end do

  ! Check positivity
  tau_negative=.false.
  do i=1,ncell
     if(tau(i)<=0.)tau_negative=.true.
  end do
  if (tau_negative) then
     write(*,*)'ERROR in solve_cooling :'
     write(*,*)'Initial temperature is negative'
     STOP
  endif

  ! Loop over active cells
  iter=0
  n=ncell
  do while(n>0)

     iter=iter+1
     if (iter > 500) then
        write(*,*) 'Too many iterations in solve_cooling',iter,n
        do i=1,n
           write(*,*)i,tau(ind(i)),T2(ind(i)),nH(ind(i)),i_nH(ind(i))
        end do
        STOP
     endif

     n_active=0
     do i=1,n
        facT=log10(tau(ind(i)))

        if(facT.le.logT2max)then

           i_T2=MIN(MAX(int((facT-c%table%T2(1))*dlog_T2)+1,1),c%table%n2-1)
           yy=facT-c%table%T2(i_T2)
           yy2=yy*yy
           yy3=yy2*yy

           ! Cooling
           fa=c%table%cool(i_nH(ind(i)),i_T2  )*w1H(ind(i))+c%table%cool(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fb=c%table%cool(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+c%table%cool(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           fprimea=c%table%cool_prime(i_nH(ind(i)),i_T2  )*w1H(ind(i))+c%table%cool_prime(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fprimeb=c%table%cool_prime(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+c%table%cool_prime(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           cool=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           cool_prime=cool/tau(ind(i))*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Heating
           fa=c%table%heat(i_nH(ind(i)),i_T2  )*w1H(ind(i))+c%table%heat(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fb=c%table%heat(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+c%table%heat(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           fprimea=c%table%heat_prime(i_nH(ind(i)),i_T2  )*w1H(ind(i))+c%table%heat_prime(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fprimeb=c%table%heat_prime(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+c%table%heat_prime(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           heat=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           heat_prime=heat/tau(ind(i))*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Compton cooling
           fa=c%table%cool_com(i_nH(ind(i)),i_T2  )*w1H(ind(i))+c%table%cool_com(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fb=c%table%cool_com(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+c%table%cool_com(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           fprimea=c%table%cool_com_prime(i_nH(ind(i)),i_T2  )*w1H(ind(i))+c%table%cool_com_prime(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fprimeb=c%table%cool_com_prime(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+c%table%cool_com_prime(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           cool_com=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           cool_com_prime=cool_com/tau(ind(i))*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Compton heating
           fa=c%table%heat_com(i_nH(ind(i)),i_T2  )*w1H(ind(i))+c%table%heat_com(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fb=c%table%heat_com(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+c%table%heat_com(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           fprimea=c%table%heat_com_prime(i_nH(ind(i)),i_T2  )*w1H(ind(i))+c%table%heat_com_prime(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fprimeb=c%table%heat_com_prime(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+c%table%heat_com_prime(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           heat_com=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           heat_com_prime=heat_com/tau(ind(i))*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Metal cooling
           fa=c%table%metal(i_nH(ind(i)),i_T2  )*w1H(ind(i))+c%table%metal(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fb=c%table%metal(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+c%table%metal(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           fprimea=c%table%metal_prime(i_nH(ind(i)),i_T2  )*w1H(ind(i))+c%table%metal_prime(i_nH(ind(i))+1,i_T2  )*w2H(ind(i))
           fprimeb=c%table%metal_prime(i_nH(ind(i)),i_T2+1)*w1H(ind(i))+c%table%metal_prime(i_nH(ind(i))+1,i_T2+1)*w2H(ind(i))
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           metal=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           metal_prime=metal/tau(ind(i))*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Total net cooling
           lambda=cool+zzz(ind(i))*metal-heat+(cool_com-heat_com)/nH(ind(i))
           lambda_prime=cool_prime+zzz(ind(i))*metal_prime-heat_prime+(cool_com_prime-heat_com_prime)/nH(ind(i))

        else

           lambda=1.42d0*1d-27*sqrt(tau(ind(i)))*1.1d0
           lambda_prime=lambda/2d0/tau(ind(i))

        endif

        wcool=MAX(abs(lambda)/tau(ind(i))*varmax,wmax(ind(i)),-lambda_prime*varmax)

        tau_old(ind(i))=tau(ind(i))
        tau(ind(i))=tau(ind(i))*(1d0+lambda_prime/wcool-lambda/tau(ind(i))/wcool)/(1d0+lambda_prime/wcool)
        time_old(ind(i))=time(ind(i))
        time(ind(i))=time(ind(i))+1d0/wcool

!!$        if(i==1)then
!!$           write(10,'(I5,10(1PE10.3,1X))')iter,tau_old(ind(i)),cool+zzz(ind(i))*metal,heat,lambda
!!$        endif

        if(time(ind(i))<time_max(ind(i)))then
           n_active=n_active+1
           ind(n_active)=ind(i)
        end if

     end do
     n=n_active
  end do
  ! End loop over active cells

  ! Compute exact time solution
  do i=1,ncell
     tau(i)=tau(i)*(time_max(i)-time_old(i))/(time(i)-time_old(i))+tau_old(i)*(time(i)-time_max(i))/(time(i)-time_old(i))
  end do

  ! Check positivity
  tau_negative=.false.
  do i=1,ncell
     if (tau(i)<=0.)tau_negative=.true.
  end do
  if (tau_negative) then
     write(*,*)'ERROR in solve_cooling :'
     write(*,*)'Final temperature is negative'
     STOP
  endif

  ! Compute delta T
  do i=1,ncell
     deltaT2(i)=tau(i)-tau_ini(i)
  end do

end subroutine solve_cooling
!=======================================================================
function J0simple(c,aexp)
!=======================================================================
! Le J0 dans le cas teyssier ou theuns
!=======================================================================
  type(cooling_t)::c
  real(kind=8) :: J0simple,aexp
  if (aexp .lt. 1d0/(1d0+c%zreioniz)) then
     J0simple=0d0
  elseif (aexp .lt. 1d0/4d0)then
     J0simple=4d0*aexp
  elseif (aexp .lt. 1d0/3d0)then
     J0simple=1d0
  else
     J0simple=1d0/(3d0*aexp)**3
  endif
  if(c%force_j0_one)J0simple=1.0d0
  J0simple=max(J0simple*c%J0in,c%J0min)
  return
end function J0simple
!=======================================================================
subroutine cmp_table(c,nH_min,nH_max,T2_min,T2_max,nbin_n,nbin_T,aexp)
!=======================================================================
  implicit none
  type(cooling_t)::c
  real(kind=8)::nH_min,nH_max,T2_min,T2_max,aexp,tmp
  integer::nbin_n,nbin_T
  integer::i_n,i_T
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec

  if(.not.c%first)then
     deallocate(c%table%cool)
     deallocate(c%table%heat)
     deallocate(c%table%cool_com)
     deallocate(c%table%heat_com)
     deallocate(c%table%metal)
     deallocate(c%table%cool_prime)
     deallocate(c%table%heat_prime)
     deallocate(c%table%cool_com_prime)
     deallocate(c%table%heat_com_prime)
     deallocate(c%table%metal_prime)
     deallocate(c%table%mu)
     deallocate(c%table%T2)
     deallocate(c%table%nH)
     if (if_species_abundances) deallocate(c%table%n_spec)
  else
     c%first=.false.
  endif

  c%table%n1=nbin_n
  c%table%n2=nbin_T
  allocate(c%table%cool(nbin_n,nbin_T))
  allocate(c%table%heat(nbin_n,nbin_T))
  allocate(c%table%cool_com(nbin_n,nbin_T))
  allocate(c%table%heat_com(nbin_n,nbin_T))
  allocate(c%table%metal(nbin_n,nbin_T))
  allocate(c%table%cool_prime(nbin_n,nbin_T))
  allocate(c%table%heat_prime(nbin_n,nbin_T))
  allocate(c%table%cool_com_prime(nbin_n,nbin_T))
  allocate(c%table%heat_com_prime(nbin_n,nbin_T))
  allocate(c%table%metal_prime(nbin_n,nbin_T))
  allocate(c%table%mu  (nbin_n,nbin_T))
  allocate(c%table%nH  (1:nbin_n))
  allocate(c%table%T2  (1:nbin_T))
  if (if_species_abundances) allocate(c%table%n_spec(nbin_n,nbin_T,1:6))

  do i_n=1,nbin_n
     tmp=log10(nH_min)+(dble(i_n)-1d0)/(dble(nbin_n)-1d0)*(log10(nH_max)-log10(nH_min))
     c%table%nH(i_n)=tmp
  end do
  do i_T=1,nbin_T
     tmp=log10(T2_min)+(dble(i_T)-1d0)/(dble(nbin_T)-1d0)*(log10(T2_max)-log10(T2_min))
     c%table%T2(i_T)=tmp
  end do

  ! Compute radiative ionization and heating rates
  call set_rates(c,t_rad_spec,h_rad_spec,aexp)

  ! Record the reference rate (HI photo-ionization rate)
  c%uvb_last_update = t_rad_spec(HI)

  ! Create the c%table
  c%table%mu=0
  c%table%cool=0
  c%table%heat=0
  c%table%cool_com=0
  c%table%heat_com=0
  c%table%metal=0
  c%table%cool_prime=0
  c%table%heat_prime=0
  c%table%cool_com_prime=0
  c%table%heat_com_prime=0
  c%table%metal_prime=0
  if (if_species_abundances) c%table%n_spec=0
  do i_n = 1,nbin_n
     call iterate(c,i_n,t_rad_spec,h_rad_spec,nbin_T,aexp)
  end do

end subroutine cmp_table
!=======================================================================
subroutine set_rates(c,t_rad_spec,h_rad_spec,aexp)
!=======================================================================
  implicit none
  type(cooling_t)::c
  real(kind=8),dimension(1:3) :: t_rad_spec,h_rad_spec
  real(kind=8) :: J0,z,aexp

  z=1d0/aexp-1d0
  if (c%high_z_realistic_ne) c%J0min=c%J0min_ref/(aexp/c%aexp_ref)**2
  if (c%teyssier) then
     J0=J0simple(c,aexp)
     t_rad_spec(HI  ) = taux_rad(c,HI  ,J0)
     t_rad_spec(HEI ) = taux_rad(c,HEI ,J0)
     t_rad_spec(HEII) = taux_rad(c,HEII,J0)
     h_rad_spec(HI  ) = heat_rad(c,HI  ,J0)
     h_rad_spec(HEI ) = heat_rad(c,HEI ,J0)
     h_rad_spec(HEII) = heat_rad(c,HEII,J0)
  elseif (c%theuns) then
     J0=J0simple(c,aexp)
     t_rad_spec(HI  ) = taux_rad_theuns(c,HI  ,J0)
     t_rad_spec(HEI ) = taux_rad_theuns(c,HEI ,J0)
     t_rad_spec(HEII) = taux_rad_theuns(c,HEII,J0)
     h_rad_spec(HI  ) = heat_rad_theuns(c,HI  ,J0)
     h_rad_spec(HEI ) = heat_rad_theuns(c,HEI ,J0)
     h_rad_spec(HEII) = heat_rad_theuns(c,HEII,J0)
  elseif (c%madau) then
     z=1d0/aexp-1d0
     t_rad_spec(HI  ) = taux_rad_madau(c,HI  ,z)
     t_rad_spec(HEI ) = taux_rad_madau(c,HEI ,z)
     t_rad_spec(HEII) = taux_rad_madau(c,HEII,z)
     h_rad_spec(HI  ) = heat_rad_madau(c,HI  ,z)
     h_rad_spec(HEI ) = heat_rad_madau(c,HEI ,z)
     h_rad_spec(HEII) = heat_rad_madau(c,HEII,z)
  elseif (c%courty) then
     t_rad_spec(HI  ) = taux_rad_courty(c,HI  ,z)
     t_rad_spec(HEI ) = taux_rad_courty(c,HEI ,z)
     t_rad_spec(HEII) = taux_rad_courty(c,HEII,z)
     h_rad_spec(HI  ) = heat_rad_courty(c,HI  ,z)
     h_rad_spec(HEI ) = heat_rad_courty(c,HEI ,z)
     h_rad_spec(HEII) = heat_rad_courty(c,HEII,z)
  endif
end subroutine set_rates
!=======================================================================
subroutine iterate(c,i_n,t_rad_spec,h_rad_spec,nbin_T,aexp)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer :: i_n
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec
  real(kind=8) :: aexp
  integer::nbin_T
  integer::i_T
  real(kind=8) ::T2,T2_eps,nH
  real(kind=8) ::mu,mu_eps
  real(kind=8) ::cool_tot,heat_tot,cool_com,heat_com,metal_tot,metal_prime
  real(kind=8) ::cool_tot_eps,heat_tot_eps,cool_com_eps,heat_com_eps
  real(kind=8),dimension(1:6) :: n_spec,n_spec_eps

  nH=10d0**c%table%nH(i_n)
  do i_T = 1,nbin_T
     T2=10d0**c%table%T2(i_T)
     ! Compute cooling, heating and mean molecular weight
     call cmp_cooling(c,T2,nH,t_rad_spec,h_rad_spec,cool_tot,heat_tot,cool_com,heat_com,mu,aexp,n_spec)
     c%table%cool(i_n,i_T)=log10(cool_tot)
     c%table%heat(i_n,i_T)=log10(heat_tot)
     c%table%cool_com(i_n,i_T)=log10(cool_com)
     c%table%heat_com(i_n,i_T)=log10(heat_com)
     c%table%mu(i_n,i_T)=mu
     if (if_species_abundances)then
        c%table%n_spec(i_n,i_T,1:6)=log10(n_spec(1:6))
     endif
     ! Compute cooling and heating derivatives
     T2_eps=10d0**(c%table%T2(i_T)+0.01d0)
     call cmp_cooling(c,T2_eps,nH,t_rad_spec,h_rad_spec,cool_tot_eps,heat_tot_eps,cool_com_eps,heat_com_eps,mu_eps,aexp,n_spec_eps)
     c%table%cool_prime(i_n,i_T)=(log10(cool_tot_eps)-log10(cool_tot))/0.01d0
     c%table%heat_prime(i_n,i_T)=(log10(heat_tot_eps)-log10(heat_tot))/0.01d0
     c%table%cool_com_prime(i_n,i_T)=(log10(cool_com_eps)-log10(cool_com))/0.01d0
     c%table%heat_com_prime(i_n,i_T)=(log10(heat_com_eps)-log10(heat_com))/0.01d0
     ! Compute metal contribution for solar metallicity
     call cmp_metals(c,T2,nH,mu,metal_tot,metal_prime,aexp)
     c%table%metal(i_n,i_T)=log10(metal_tot)
     c%table%metal_prime(i_n,i_T)=metal_prime
  end do
end subroutine iterate
!=======================================================================
subroutine cmp_metals(c,T2,nH,mu,metal_tot,metal_prime,aexp)
!=======================================================================
  implicit none
  type(cooling_t)::c
  real(kind=8) ::T2,nH,mu,metal_tot,metal_prime,aexp

  ! Compute cooling enhancement due to metals
  ! Cloudy at solar metalicity
  real(kind=8),dimension(1:91) :: temperature_cc07 = (/ &
       & 3.9684,4.0187,4.0690,4.1194,4.1697,4.2200,4.2703, &
       & 4.3206,4.3709,4.4212,4.4716,4.5219,4.5722,4.6225, &
       & 4.6728,4.7231,4.7734,4.8238,4.8741,4.9244,4.9747, &
       & 5.0250,5.0753,5.1256,5.1760,5.2263,5.2766,5.3269, &
       & 5.3772,5.4275,5.4778,5.5282,5.5785,5.6288,5.6791, &
       & 5.7294,5.7797,5.8300,5.8804,5.9307,5.9810,6.0313, &
       & 6.0816,6.1319,6.1822,6.2326,6.2829,6.3332,6.3835, &
       & 6.4338,6.4841,6.5345,6.5848,6.6351,6.6854,6.7357, &
       & 6.7860,6.8363,6.8867,6.9370,6.9873,7.0376,7.0879, &
       & 7.1382,7.1885,7.2388,7.2892,7.3395,7.3898,7.4401, &
       & 7.4904,7.5407,7.5911,7.6414,7.6917,7.7420,7.7923, &
       & 7.8426,7.8929,7.9433,7.9936,8.0439,8.0942,8.1445, &
       & 8.1948,8.2451,8.2955,8.3458,8.3961,8.4464,8.4967 /)
! Cooling from metals only (without the contribution of H and He)
! log cooling rate in [erg s-1 cm3]
! S. Ploeckinger 06/2015
  real(kind=8),dimension(1:91) :: excess_cooling_cc07 = (/ &
       &  -24.9082, -24.9082, -24.5503, -24.0898, -23.5328, -23.0696, -22.7758, &
       &  -22.6175, -22.5266, -22.4379, -22.3371, -22.2289, -22.1181, -22.0078, &
       &  -21.8992, -21.7937, -21.6921, -21.5961, -21.5089, -21.4343, -21.3765, &
       &  -21.3431, -21.3274, -21.3205, -21.3142, -21.3040, -21.2900, -21.2773, &
       &  -21.2791, -21.3181, -21.4006, -21.5045, -21.6059, -21.6676, -21.6877, &
       &  -21.6934, -21.7089, -21.7307, -21.7511, -21.7618, -21.7572, -21.7532, &
       &  -21.7668, -21.7860, -21.8129, -21.8497, -21.9035, -21.9697, -22.0497, &
       &  -22.1327, -22.2220, -22.3057, -22.3850, -22.4467, -22.4939, -22.5205, &
       &  -22.5358, -22.5391, -22.5408, -22.5408, -22.5475, -22.5589, -22.5813, &
       &  -22.6122, -22.6576, -22.7137, -22.7838, -22.8583, -22.9348, -23.0006, &
       &  -23.0547, -23.0886, -23.1101, -23.1139, -23.1147, -23.1048, -23.1017, &
       &  -23.0928, -23.0969, -23.0968, -23.1105, -23.1191, -23.1388, -23.1517, &
       &  -23.1717, -23.1837, -23.1986, -23.2058, -23.2134, -23.2139, -23.2107 /)
  real(kind=8),dimension(1:91) :: excess_prime_cc07 = (/ &
       &   2.0037,  4.7267, 12.2283, 13.5820,  9.8755,  4.8379,  1.8046, &
       &   1.4574,  1.8086,  2.0685,  2.2012,  2.2250,  2.2060,  2.1605, &
       &   2.1121,  2.0335,  1.9254,  1.7861,  1.5357,  1.1784,  0.7628, &
       &   0.1500, -0.1401,  0.1272,  0.3884,  0.2761,  0.1707,  0.2279, &
       &  -0.2417, -1.7802, -3.0381, -2.3511, -0.9864, -0.0989,  0.1854, &
       &  -0.1282, -0.8028, -0.7363, -0.0093,  0.3132,  0.1894, -0.1526, &
       &  -0.3663, -0.3873, -0.3993, -0.6790, -1.0615, -1.4633, -1.5687, &
       &  -1.7183, -1.7313, -1.8324, -1.5909, -1.3199, -0.8634, -0.5542, &
       &  -0.1961, -0.0552,  0.0646, -0.0109, -0.0662, -0.2539, -0.3869, &
       &  -0.6379, -0.8404, -1.1662, -1.3930, -1.6136, -1.5706, -1.4266, &
       &  -1.0460, -0.7244, -0.3006, -0.1300,  0.1491,  0.0972,  0.2463, &
       &   0.0252,  0.1079, -0.1893, -0.1033, -0.3547, -0.2393, -0.4280, &
       &  -0.2735, -0.3670, -0.2033, -0.2261, -0.0821, -0.0754,  0.0634 /)
  real(kind=8),dimension(1:50)::z_courty=(/ &
       & 0.00000,0.04912,0.10060,0.15470,0.21140,0.27090,0.33330,0.39880, &
       & 0.46750,0.53960,0.61520,0.69450,0.77780,0.86510,0.95670,1.05300, &
       & 1.15400,1.25900,1.37000,1.48700,1.60900,1.73700,1.87100,2.01300, &
       & 2.16000,2.31600,2.47900,2.64900,2.82900,3.01700,3.21400,3.42100, &
       & 3.63800,3.86600,4.10500,4.35600,4.61900,4.89500,5.18400,5.48800, &
       & 5.80700,6.14100,6.49200,6.85900,7.24600,7.65000,8.07500,8.52100, &
       & 8.98900,9.50000 /)
  real(kind=8),dimension(1:50)::phi_courty=(/ &
       & 0.0499886,0.0582622,0.0678333,0.0788739,0.0915889,0.1061913,0.1229119, &
       & 0.1419961,0.1637082,0.1883230,0.2161014,0.2473183,0.2822266,0.3210551, &
       & 0.3639784,0.4111301,0.4623273,0.5172858,0.5752659,0.6351540,0.6950232, &
       & 0.7529284,0.8063160,0.8520859,0.8920522,0.9305764,0.9682031,1.0058810, &
       & 1.0444020,1.0848160,1.1282190,1.1745120,1.2226670,1.2723200,1.3231350, &
       & 1.3743020,1.4247480,1.4730590,1.5174060,1.5552610,1.5833640,1.5976390, &
       & 1.5925270,1.5613110,1.4949610,1.3813710,1.2041510,0.9403100,0.5555344, &
       & 0.0000000 /)
  real(kind=8)::TT,lTT,deltaT,lcool1,lcool2,lcool1_prime,lcool2_prime
  real(kind=8)::ZZ,deltaZ
  real(kind=8)::c1=0.4d0,c2=10.0d0,TT0=1d5,TTC=1d6,alpha1=0.15d0
  real(kind=8)::ux,g_courty,f_courty=1d0,g_courty_prime,f_courty_prime
  integer::iT,iZ

  ZZ=1d0/aexp-1d0
  TT=T2*mu
  lTT=log10(TT)

  ! This is a simple model to take into account the ionization background
  ! on metal cooling (calibrated using CLOUDY).
  if(c%madau.or.c%courty)then
     if(ZZ.le.0.0.or.ZZ.ge.z_courty(50))then
        ux=0
     else
        iZ=1
        do while(ZZ>z_courty(iZ+1).and.iZ<49)
           iZ=iZ+1
        end do
        deltaZ=z_courty(iZ+1)-z_courty(iZ)
        ux=1d-4*(phi_courty(iZ+1)*(ZZ-z_courty(iZ))/deltaZ &
             & + phi_courty(iZ)*(z_courty(iZ+1)-ZZ)/deltaZ )/nH
     endif
  else ! Theuns or Teyssier
     ux=1d-4*J0simple(c,aexp)/1d-22/nH
  endif
  g_courty=c1*(TT/TT0)**alpha1+c2*exp(-TTC/TT)
  g_courty_prime=(c1*alpha1*(TT/TT0)**alpha1+c2*exp(-TTC/TT)*TTC/TT)/TT
  f_courty=1d0/(1d0+ux/g_courty)
  f_courty_prime=ux/g_courty/(1d0+ux/g_courty)**2*g_courty_prime/g_courty

  if(lTT.ge.temperature_cc07(91))then
     metal_tot=1d-100
     metal_prime=0
  else if(lTT.ge.1.0)then
     lcool1=-100
     lcool1_prime=0
     if(lTT.ge.temperature_cc07(1))then
        iT=1+int((lTT-temperature_cc07(1))/(temperature_cc07(91)-temperature_cc07(1))*90)
        iT=min(iT,90)
        iT=max(iT,1)
        deltaT=temperature_cc07(iT+1)-temperature_cc07(iT)
        lcool1 = excess_cooling_cc07(iT+1)*(lTT-temperature_cc07(iT))/deltaT  &
             & + excess_cooling_cc07(iT)*(temperature_cc07(iT+1)-lTT)/deltaT
        lcool1_prime  = excess_prime_cc07(iT+1)*(lTT-temperature_cc07(iT))/deltaT  &
                    & + excess_prime_cc07(iT)*(temperature_cc07(iT+1)-lTT)/deltaT
     endif
     ! Fine structure cooling from infrared lines
     lcool2=-31.522879d0+2*lTT-20.0d0/TT-TT*4.342944d-5
     lcool2_prime=2d0+(20/TT-TT*4.342944d-5)*log(10d0)
     ! Total metal cooling and temperature derivative
     metal_tot=10d0**lcool1+10d0**lcool2
     metal_prime=(10d0**lcool1*lcool1_prime+10d0**lcool2*lcool2_prime)/metal_tot
     metal_prime=metal_prime*f_courty+metal_tot*f_courty_prime
     metal_tot=metal_tot*f_courty
  else
     metal_tot=1d-100
     metal_prime=0
  endif

end subroutine cmp_metals
!=======================================================================
subroutine cmp_cooling(c,T2,nH,t_rad_spec,h_rad_spec,cool_tot,heat_tot,cool_com,heat_com,mu_out,aexp,n_spec)
!=======================================================================
  implicit none
  type(cooling_t)::c
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec
  real(kind=8) ::T2,nH,cool_tot,heat_tot,cool_com,heat_com,mu_out,aexp
  real(kind=8) ::mu,mu_old,err_mu,mu_left,mu_right
  real(kind=8) ::T
  real(kind=8) ::n_E,n_HI,n_HII,n_HEI,n_HEII,n_HEIII,n_TOT
  real(kind=8),dimension(1:6)::n_spec
  real(kind=8) ::cb1,cb2,cb3,ci1,ci2,ci3,cr1,cr2,cr3,cd,ce1,ce2,ce3,ch1,ch2,ch3,coc,coh
  integer::niter

  ! Iteration to find mu
  err_mu=1
  mu_left=0.5d0
  mu_right=1.3d0
  niter=0
  do while (err_mu > 1d-4 .and. niter <= 50)
     mu_old=0.5d0*(mu_left+mu_right)
     T = T2*mu_old
     call cmp_chem_eq(c,T,nH,t_rad_spec,n_spec,n_TOT,mu)
     err_mu = (mu-mu_old)/mu_old
     if(err_mu>0.)then
        mu_left =0.5d0*(mu_left+mu_right)
        mu_right=mu_right
     else
        mu_left =mu_left
        mu_right=0.5d0*(mu_left+mu_right)
     end if
     err_mu=ABS(err_mu)
     niter=niter+1
  end do
  if (niter > 50) then
     write(*,*) 'ERROR in cmp_cooling : too many iterations.'
     STOP
  endif

  ! Get equilibrium abundances
  n_E     = n_spec(1) ! electrons
  n_HI    = n_spec(2) ! H
  n_HII   = n_spec(3) ! H+
  n_HEI   = n_spec(4) ! He
  n_HEII  = n_spec(5) ! He+
  n_HEIII = n_spec(6) ! He++
  ! Bremstrahlung
  cb1 = cool_bre(HI  ,T)*n_E*n_HII  /nH**2
  cb2 = cool_bre(HEI ,T)*n_E*n_HEII /nH**2
  cb3 = cool_bre(HEII,T)*n_E*n_HEIII/nH**2
  ! Ionization cooling
  ci1 = cool_ion(c,HI  ,T)*n_E*n_HI   /nH**2
  ci2 = cool_ion(c,HEI ,T)*n_E*n_HEI  /nH**2
  ci3 = cool_ion(c,HEII,T)*n_E*n_HEII /nH**2
  ! Recombination cooling
  cr1 = cool_rec(HI  ,T)*n_E*n_HII  /nH**2
  cr2 = cool_rec(HEI ,T)*n_E*n_HEII /nH**2
  cr3 = cool_rec(HEII,T)*n_E*n_HEIII/nH**2
  ! Dielectric recombination cooling
  cd  = cool_die(     T)*n_E*n_HEII /nH**2
  ! Line cooling
  ce1 = cool_exc(HI  ,T)*n_E*n_HI   /nH**2
  ce2 = cool_exc(HEI, T)*n_E*n_HEI  /nH**2
  ce3 = cool_exc(HEII,T)*n_E*n_HEII /nH**2
  ! Radiative heating
  ch1 = h_rad_spec(HI  )    *n_HI   /nH**2
  ch2 = h_rad_spec(HEI )    *n_HEI  /nH**2
  ch3 = h_rad_spec(HEII)    *n_HEII /nH**2
  ! Total cooling and heating rates
  heat_tot = ch1+ch2+ch3
  cool_tot = cb1+cb2+cb3+ci1+ci2+ci3+cr1+cr2+cr3+cd+ce1+ce2+ce3
  ! Compton cooling
  coc = cool_compton(T,aexp)*n_E/nH
  cool_com = coc
  ! Compton heating
  coh = heat_compton(T,aexp)*n_E/nH
  heat_com = coh
  ! Mean molecular weight
  mu_out = mu

end subroutine cmp_cooling
!=======================================================================
subroutine cmp_chem_eq(c,T,n_H,t_rad_spec,n_spec,n_TOT,mu)
!=======================================================================
  implicit none
  type(cooling_t)::c
  real(kind=8)::T,n_H,n_TOT,mu
  real(kind=8),dimension(1:3)::t_rad_spec
  real(kind=8),dimension(1:6)::n_spec
  real(kind=8)::xx,yy
  real(kind=8)::n_HI,n_HII,n_HEI,n_HEII,n_HEIII,n_E
  real(kind=8)::t_rad_HI,t_rad_HEI,t_rad_HEII
  real(kind=8)::t_rec_HI,t_rec_HEI,t_rec_HEII
  real(kind=8)::t_ion_HI,t_ion_HEI,t_ion_HEII
  real(kind=8)::t_ion2_HI,t_ion2_HEI,t_ion2_HEII
  real(kind=8)::x1,err_nE

  xx=(1d0-c%Y)
  yy=c%Y/(1d0-c%Y)/4d0

  t_rad_HI   = t_rad_spec(HI)
  t_rec_HI   = taux_rec  (c,HI,T)
  t_ion_HI   = taux_ion  (c,HI,T)

  t_rad_HEI  = t_rad_spec(HEI)
  t_rec_HEI  = taux_rec  (c,HEI,T)
  t_ion_HEI  = taux_ion  (c,HEI,T)

  t_rad_HEII = t_rad_spec(HEII)
  t_rec_HEII = taux_rec  (c,HEII,T)
  t_ion_HEII = taux_ion  (c,HEII,T)

  n_E = n_H
  err_nE = 1

  do while(err_nE > 1d-8)

     t_ion2_HI   = t_ion_HI   + t_rad_HI  /MAX(n_E,1d-15*n_H)
     t_ion2_HEI  = t_ion_HEI  + t_rad_HEI /MAX(n_E,1d-15*n_H)
     t_ion2_HEII = t_ion_HEII + t_rad_HEII/MAX(n_E,1d-15*n_H)

     n_HI  = t_rec_HI/(t_ion2_HI+t_rec_HI)*n_H
     n_HII = t_ion2_HI/(t_ion2_HI+t_rec_HI)*n_H

     x1 = (t_rec_HEII*t_rec_HEI+t_ion2_HEI*t_rec_HEII+t_ion2_HEII*t_ion2_HEI)

     n_HEIII = yy*t_ion2_HEII*t_ion2_HEI/x1*n_H
     n_HEII  = yy*t_ion2_HEI *t_rec_HEII/x1*n_H
     n_HEI   = yy*t_rec_HEII *t_rec_HEI /x1*n_H

     err_nE = ABS((n_E - (n_HII + n_HEII + 2*n_HEIII))/n_H)
     n_E = 0.5d0*n_E+0.5d0*(n_HII + n_HEII + 2*n_HEIII)

  end do

  n_TOT    =n_E+n_HI+n_HII+n_HEI+n_HEII+n_HEIII
  mu       =n_H/xx/n_TOT
  n_spec(1)=n_E
  n_spec(2)=n_HI
  n_spec(3)=n_HII
  n_spec(4)=n_HEI
  n_spec(5)=n_HEII
  n_spec(6)=n_HEIII

end subroutine cmp_chem_eq
!=======================================================================
function cool_bre(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,cool_bre
  cool_bre = 0
  if(ispec==HI  )cool_bre = 1.42D-27*sqrt(T)*(1.1D0+0.34D0*exp(-(5.5D0-log10(T))**2 /3d0))
  if(ispec==HEI )cool_bre = 1.42D-27*sqrt(T)*(1.1D0+0.34D0*exp(-(5.5D0-log10(T))**2 /3d0))
  if(ispec==HEII)cool_bre = 5.68D-27*sqrt(T)*(1.1D0+0.34D0*exp(-(5.5D0-log10(T))**2 /3d0))
  return
end function cool_bre
!=======================================================================
function cool_exc(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,cool_exc,T5
  T5=1d-5*T
  cool_exc = 0
  if(ispec==HI  )cool_exc = 7.50D-19/(1d0+sqrt(T5))              *exp(-118348d0/T)
  if(ispec==HEI )cool_exc = 9.10D-27/(1d0+sqrt(T5))/(T**0.1687D0)*exp(-13179d0/T)
  if(ispec==HEII)cool_exc = 5.54D-17/(1d0+sqrt(T5))/(T**0.397D0 )*exp(-473638d0/T)
  return
end function cool_exc
!=======================================================================
function cool_rec(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,cool_rec
  real(kind=8)   ::T3, T6
  T3 = 1d-03*T
  T6 = 1d-06*T
  cool_rec = 0
  if(ispec==HI  )cool_rec = 8.70D-27*SQRT(T)/T3**(0.2d0)/(1d0+T6**0.7D0)
  if(ispec==HEI )cool_rec = 1.55D-26*T**0.3647D0
  if(ispec==HEII)cool_rec = 3.48D-26*SQRT(T)/T3**(0.2d0)/(1d0+T6**0.7D0)
  return
end function cool_rec
!=======================================================================
function cool_die(T)
!=======================================================================
  implicit none
  real(kind=8) :: T,cool_die
  cool_die=1.24D-13*T**(-1.5D0)*exp(-470000d0/T)*(1d0+0.3D0*exp(-94000d0/T))
  return
end function cool_die
!=======================================================================
function taux_rec(c,ispec,T)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer::ispec
  real(kind=8)   ::T,taux_rec
  real(kind=8)   ::T3, T6
  T3 = 1d-03*T
  T6 = 1d-06*T
  taux_rec = 0
  if(ispec==HI  )taux_rec = c%dumfac_rec*8.40d-11/SQRT(T)/T3**(0.2d0)/(1d0+T6**0.7d0)
  if(ispec==HEI )taux_rec = 1.50d-10/T**0.6353d0+taux_die(T)
  if(ispec==HEII)taux_rec = 3.36d-10/SQRT(T)/T3**(0.2d0)/(1d0+T6**0.7d0)
  return
end function taux_rec
!=======================================================================
function taux_die(T)
!=======================================================================
  implicit none
  real(kind=8) :: T,taux_die
  taux_die=1.9D-3*T**(-1.5D0)*exp(-470000d0/T)*(1d0+0.3D0*exp(-94000d0/T))
  return
end function taux_die
!=======================================================================
function cool_ion(c,ispec,T)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer::ispec
  real(kind=8)   ::T,cool_ion
  real(kind=8)   ::T5
  T5 = 1d-05*T
  cool_ion = 0
  if(ispec==HI  )cool_ion = c%dumfac_ion*1.27D-21*SQRT(T)/(1d0+SQRT(T5))*EXP(-157809.1D0/T)
  if(ispec==HEI )cool_ion = c%dumfac_ion*9.38D-22*SQRT(T)/(1d0+SQRT(T5))*EXP(-285335.4D0/T)
  if(ispec==HEII)cool_ion = c%dumfac_ion*4.95D-22*SQRT(T)/(1d0+SQRT(T5))*EXP(-631515.0D0/T)
  return
end function cool_ion
!=======================================================================
function cool_compton(T,aexp)
!=======================================================================
  implicit none
  real(kind=8) ::T,aexp,cool_compton
  cool_compton=5.406D-36*T/aexp**4
  return
end function cool_compton
!=======================================================================
function heat_compton(T,aexp)
!=======================================================================
  implicit none
  real(kind=8) ::T,aexp,heat_compton
  real(kind=8) ::T5
  T5 = 1d-05*T
  heat_compton=5.406D-36*2.726D0/aexp**5
  return
end function heat_compton
!=======================================================================
function taux_ion(c,ispec,T)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer::ispec
  real(kind=8)   :: T,taux_ion
  real(kind=8)   :: T5
  T5 = 1d-05*T
  taux_ion = 0
  if(ispec==HI  )taux_ion = c%dumfac_ion*5.85D-11*SQRT(T)/(1d0+SQRT(T5))*EXP(-157809.1D0/T)
  if(ispec==HEI )taux_ion = c%dumfac_ion*2.38D-11*SQRT(T)/(1d0+SQRT(T5))*EXP(-285335.4D0/T)
  if(ispec==HEII)taux_ion = c%dumfac_ion*5.68D-12*SQRT(T)/(1d0+SQRT(T5))*EXP(-631515.0D0/T)
  return
end function taux_ion
!=======================================================================
function J_nu(c,e,J0)
!=======================================================================
  implicit none
  type(cooling_t)::c
  real(kind=8) :: e,J_nu,e_L,J0,Jloc
  Jloc = max(J0,c%J0min)
  e_L  = 13.598d0*eV2erg
  J_nu = Jloc*(e_L/e)
  return
end function J_nu
!=======================================================================
function sigma_rad(e,ispec)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)::sigma_rad,e,e_i,xxx,alph
  if(ispec==HI  )e_i = 13.598D0*eV2erg
  if(ispec==HEI )e_i = 24.587D0*eV2erg
  if(ispec==HEII)e_i = 54.416D0*eV2erg
  xxx = e/e_i
  alph = sqrt(xxx-1.0d0)
  sigma_rad=0
  if(ispec==HI  )sigma_rad = 6.30D-18/xxx**4*exp(4d0-4d0*atan(alph)/alph) &
       &                             /(1d0-exp(-twopi/alph))
  if(ispec==HEI )sigma_rad = 7.42D-18*(1.66D0/xxx**2.05D0-0.66D0/xxx**3.05D0)
  if(ispec==HEII)sigma_rad = 1.58D-18/xxx**4*exp(4d0-4d0*atan(alph)/alph) &
       &                             /(1d0-exp(-twopi/alph))
  return
end function sigma_rad
!=======================================================================
function taux_rad(c,ispec,J0)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer::ispec
  real(kind=8) :: J0,taux_rad,e_i,e,de,error,integ
  if(ispec==HI  )e_i = 13.598D0*eV2erg
  if(ispec==HEI )e_i = 24.587D0*eV2erg
  if(ispec==HEII)e_i = 54.416D0*eV2erg
  integ = 0.0d0
  e = e_i
  de = e/100d0
  error = 1d0
  do while(error>1d-6)
     e = e + de
     de = e/100d0
     error = 2.0d0*twopi*J_nu(c,e,J0)*sigma_rad(e,ispec)*de/e
     integ = integ + error
     error = error/abs(integ)
  end do
  taux_rad = integ/hplanck
  return
end function taux_rad
!=======================================================================
function taux_rad_madau(c,ispec,z)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer :: ispec
  real(kind=8) :: z,taux_rad_madau,tt
  taux_rad_madau=0d0
  if (z < 15d0) then
     if (ispec==HI  ) taux_rad_madau=c%normfacJ0*exp(-31.04D0+2.795D0*z-0.5589D0*z**2)
     if (ispec==HEI ) taux_rad_madau=c%normfacJ0*exp(-31.08D0+2.822D0*z-0.5664D0*z**2)
     if (ispec==HEII) taux_rad_madau=c%normfacJ0*exp(-34.30D0+1.826D0*z-0.3899D0*z**2)
  endif
  tt=taux_rad_theuns(c,ispec,c%J0min)
  taux_rad_madau=max(tt,taux_rad_madau)
  return
end function taux_rad_madau
!=======================================================================
function taux_rad_theuns(c,ispec,J0)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer :: ispec
  real(kind=8) :: J0,taux_rad_theuns
  taux_rad_theuns=0
  if (ispec==HI  ) taux_rad_theuns=1.26D10*J0/(3d0+c%alpha)
  if (ispec==HEI ) taux_rad_theuns=1.48D10*J0*0.553D0**c%alpha &
                     & *(1.66D0/(c%alpha+2.05D0)-0.66D0/(c%alpha+3.05D0))
  if (ispec==HEII) taux_rad_theuns=3.34D9*J0*0.249D0**c%alpha/(3d0+c%alpha)
  return
end function taux_rad_theuns
!=======================================================================
function taux_rad_courty(c,ispec,z)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer :: ispec,i,iweinb
  real(kind=8) :: z,zz,taux_rad_courty,hh,tt
  taux_rad_courty=0d0
  if (z < c%zreioniz) then
     if (ispec==HI  ) iweinb=1
     if (ispec==HEI ) iweinb=2
     if (ispec==HEII) iweinb=3
     hh=0d0
     zz=max(z,1.0d-15)
     do i=0,Nordercourty
        hh=hh+c%coefcourty(i,iweinb)*zz**i
     enddo
     taux_rad_courty=10**hh
  endif
  tt=taux_rad_theuns(c,ispec,c%J0min)
  taux_rad_courty=max(tt,taux_rad_courty)
  return
end function taux_rad_courty
!=======================================================================
function heat_rad(c,ispec,J0)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer::ispec
  real(kind=8) :: J0,heat_rad,e_i,e,de,error,integ
  if(ispec==HI  )e_i = 13.598D0*eV2erg
  if(ispec==HEI )e_i = 24.587D0*eV2erg
  if(ispec==HEII)e_i = 54.416D0*eV2erg
  integ = 0.0d0
  e = e_i
  de = e/100d0
  error = 1d0
  do while(error>1d-6)
     e = e + de
     de = e/100d0
     error = 2.0d0*twopi*J_nu(c,e,J0)*sigma_rad(e,ispec)*(e/e_i-1d0)*de/e
     integ = integ + error
     error=error/abs(integ)
  end do
  heat_rad = integ/hplanck*e_i
  return
end function heat_rad
!=======================================================================
function heat_rad_madau(c,ispec,z)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer :: ispec
  real(kind=8) :: z,heat_rad_madau,tt
  heat_rad_madau=0d0
  if (z < 15d0) then
     if (ispec==HI  ) heat_rad_madau=c%normfacJ0*exp(-56.62D0+2.788D0*z-0.5594D0*z**2)
     if (ispec==HEI ) heat_rad_madau=c%normfacJ0*exp(-56.06D0+2.800D0*z-0.5532D0*z**2)
     if (ispec==HEII) heat_rad_madau=c%normfacJ0*exp(-58.67D0+1.888D0*z-0.3947D0*z**2)
  endif
  tt=heat_rad_theuns(c,ispec,c%J0min)
  heat_rad_madau=max(tt,heat_rad_madau)
  return
end function heat_rad_madau
!=======================================================================
function heat_rad_theuns(c,ispec,J0)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer :: ispec
  real(kind=8) :: J0,heat_rad_theuns
  heat_rad_theuns=0
  if (ispec==HI  ) heat_rad_theuns=(2.91D-1*J0/(2d0+c%alpha))/(3d0+c%alpha)
  if (ispec==HEI ) heat_rad_theuns=5.84D-1*J0*0.553D0**c%alpha* &
                 & (1.66D0/(c%alpha+1.05D0)-2.32D0/(c%alpha+2.05D0)+0.66D0/(c%alpha+3.05D0))
  if (ispec==HEII) heat_rad_theuns=(2.92D-1*J0*0.249D0**c%alpha/(2d0+c%alpha))/(3d0+c%alpha)
  return
end function heat_rad_theuns
!=======================================================================
function heat_rad_courty(c,ispec,z)
!=======================================================================
  implicit none
  type(cooling_t)::c
  integer :: ispec,i,iweinb
  real(kind=8) :: z,zz,heat_rad_courty,hh,tt
  heat_rad_courty=0d0
  if (z < c%zreioniz) then
     if (ispec==HI  ) iweinb=4
     if (ispec==HEI ) iweinb=5
     if (ispec==HEII) iweinb=6
     hh=0d0
     zz=max(z,1.0d-15)
     do i=0,Nordercourty
        hh=hh+c%coefcourty(i,iweinb)*zz**i
     enddo
     heat_rad_courty=10**hh
  endif
  tt=heat_rad_theuns(c,ispec,c%J0min)
  heat_rad_courty=max(tt,heat_rad_courty)
  return
end function heat_rad_courty
!=======================================================================
function HsurH0(z,omega0,omegaL,OmegaR)
!=======================================================================
  implicit none
  real(kind=8) :: HsurH0,z,omega0,omegaL,omegaR
  HsurH0=sqrt(Omega0*(1d0+z)**3+OmegaR*(1d0+z)**2+OmegaL)
end function HsurH0

!=======================================================================
function cooling_needs_update(c, aexp)
!=======================================================================
  implicit none
  logical :: cooling_needs_update
  type(cooling_t), intent(in) :: c
  real(kind=8), intent(in) :: aexp
  real(kind=8), dimension(1:3) :: t_rad_spec, h_rad_spec

  ! Evaluate current rates at scale factor aexp
  call set_rates(c, t_rad_spec, h_rad_spec, aexp)

  ! Determine if an update is needed
  if (c%uvb_last_update < 0.0d0) then
     cooling_needs_update = .true.
  else
     cooling_needs_update = (abs(t_rad_spec(HI) - c%uvb_last_update) > c%uvb_delta * c%uvb_last_update)
  endif
end function cooling_needs_update

end module cooling_module
