subroutine ecophys_para(self,parsout)
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(29),intent(out) :: parsout
  real(rk) :: mumax_incr2,b_Mumax_large2,b_Mumax_small2,a_Mumax_large2,a_Mumax_small2
  real(rk) :: b_Qmin_N2,a_Qmin_N2,b_Qmax_N2,a_Qmax_N2,b_Vmax_N2,a_Vmax_N2
  real(rk) :: a_affin_N2,a_affin_P2,a_carbon2
  real(rk) :: b_affin_N2,b_affin_P2,b_carbon2
  real(rk) :: b_Qmin_P2,a_Qmin_P2,b_Qmax_P2,a_Qmax_P2,b_Vmax_P2,a_Vmax_P2

  !! Phytoplankton eco-physiological parameters
  b_qmin_N2=10.0**self%b_qmin_N
  b_qmax_N2=10.0**self%b_qmax_N
  b_vmax_N2=10.0**self%b_vmax_N
  !b_kn_N2=10.0**self%b_kn_N

  b_carbon2=10.0**self%b_carbon
  b_qmax_P2=10.0**self%b_qmax_P
  b_qmin_P2=10.0**self%b_qmin_P
  b_vmax_P2=10.0**self%b_vmax_P
  !b_kn_P2=10.0**self%b_kn_P

  b_affin_N2=10.0**self%b_affin_N
  b_affin_P2=10.0**self%b_affin_P

  !! Nonlinear mumax
  b_mumax_large2=10.0**self%b_mumax_large
  b_Mumax_large2=b_mumax_large2*(pi_over_six)**self%a_mumax_large
  a_Mumax_large2=3*self%a_mumax_large

  b_Mumax_small2=10.0**self%b_mumax_small
  b_Mumax_small2=b_mumax_small2*(pi_over_six)**self%a_mumax_small
  a_Mumax_small2=3*self%a_mumax_small

  !! Conversion of Phytoplankton eco-physiological parameters to ESD and mole-C base --- Nitorgen
  call convert_BGCparams(b_qmin_N2,self%a_qmin_N,self%a_carbon,b_carbon2,b_qmin_N2,a_qmin_N2)
  call convert_BGCparams(b_qmax_N2,self%a_qmax_N,self%a_carbon,b_carbon2,b_qmax_N2,a_qmax_N2)
  call convert_BGCparams(b_vmax_N2,self%a_vmax_N,self%a_carbon,b_carbon2,b_vmax_N2,a_vmax_N2)

  !b_Kn_N2 = b_kn_N2*(acos(-1.0)/6.)**self%a_kn_N
  !a_Kn_N2 = 3*(self%a_kn_N)

  !! Nutrient affinity, m^3 mmol-C d^1
  !a_affin_N2= a_Vmax_N2/a_Kn_N2 !-1
  !b_affin_N2= b_Vmax_N2/b_Kn_N2 !0.4

  !! Conversion of Phytoplankton eco-physiological parameters to ESD and mole-C base --- Phosphorous
  call convert_BGCparams(b_qmin_P2,self%a_qmin_P,self%a_carbon,b_carbon2,b_qmin_P2,a_qmin_P2)
  call convert_BGCparams(b_qmax_P2,self%a_qmax_P,self%a_carbon,b_carbon2,b_qmax_P2,a_qmax_P2)
  call convert_BGCparams(b_vmax_P2,self%a_vmax_P,self%a_carbon,b_carbon2,b_vmax_P2,a_vmax_P2)

  a_Qmax_P2 = 0.0_rk
  !b_Kn_P2=b_kn_P2*(acos(-1.0)/6.)**self%a_kn_P
  !a_Kn_P2=3*(self%a_kn_P)

  call convert_BGCparams(b_affin_N2,self%a_affin_N,self%a_carbon,b_carbon2,b_affin_N2,a_affin_N2)
  call convert_BGCparams(b_affin_P2,self%a_affin_P,self%a_carbon,b_carbon2,b_affin_P2,a_affin_P2)

  parsout=0.0_rk
  parsout(1) = b_Mumax_small2
  parsout(2) = a_Mumax_small2
  parsout(3) = self%mumax_incr
  parsout(4) = b_Mumax_large2
  parsout(5) = a_Mumax_large2
  parsout(6) = b_Qmin_N2
  parsout(7) = a_Qmin_N2
  parsout(8) = b_Qmax_N2
  parsout(9) = a_Qmax_N2
  parsout(10) = b_Vmax_N2
  parsout(11) = a_Vmax_N2
  parsout(12) = b_affin_N2
  parsout(13) = a_affin_N2
  !parsout(14) = b_Kn_N2
  !parsout(15) = a_Kn_N2
  parsout(16) = b_Qmin_P2
  parsout(17) = a_Qmin_P2
  parsout(18) = b_Qmax_P2
  parsout(19) = a_Qmax_P2
  parsout(20) = b_Vmax_P2
  parsout(21) = a_Vmax_P2
  !parsout(22) = b_Kn_P2
  !parsout(23) = a_Kn_P2
  parsout(24) = b_affin_P2
  parsout(25) = a_affin_P2
  !parsout(26) = b_mumax2
  !parsout(27) = a_mumax2
  parsout(28) = b_carbon2
  parsout(29) = a_carbon2
  return
end subroutine ecophys_para

!--------------------------------------------------------------------
subroutine f_T(self,meantemp,f_T_out)
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(2),intent(out) :: f_T_out

  real(rk)::ft_Phy,ft_zoo
  real(rk) :: meantemp

  fT_Phy = self%Tcons_phy**((meantemp-self%T_ref)/10._rk)
  fT_zoo = self%Tcons_zoo**((meantemp-self%T_ref)/10._rk)
  f_T_out = (/fT_Phy,fT_zoo/)

  return
end subroutine f_T


!-----------------------------------------------------------------------------
subroutine  F_Co2sr(self,pCO2,f_co2)
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),intent(in) :: pCO2
  real(rk),dimension(self%phyto_num),intent(out) :: f_co2

  f_co2 = (1.0_rk-exp(-self%a_co2*pco2))/(1.0_rk+self%a_star*exp(self%log_ESD-self%a_co2*pco2))

  return
end subroutine F_Co2sr

!-----------------------------------------------------------------------------
subroutine f_parsr(self,Phy,Q_N,f_co2,F_T,par,mixl,f_par)
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in):: Phy! Phytoplankton biomass concentration, mmol-C m^-3
  real(rk),dimension(self%phyto_num),intent(in):: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
  real(rk),dimension(self%phyto_num),intent(in):: f_co2! CO2 forcing
  real(rk),dimension(2),intent(in):: F_T! Temperature dependency for phytoplankton
  real(rk),intent (in)::par
  real(rk),intent(in) :: mixl
  real(rk),dimension(self%phyto_num),intent(out)::f_par! PAR forcing

  real(rk)::k, par_w

  k = self%kbg + sum(Q_N*Phy)*self%k_phyN
  par_w = par/(mixl*k)*(1.0_rk-exp(-1.0_rk*k*mixl)) ! par_w: Average light intensity within mixed layer depth,
  f_par = 1.0_rk-exp(-(self%a_par*par_w*Q_N)/(self%mumax*f_co2*F_T(1))) !OG

  return
end subroutine f_parsr


!-----------------------------------------------------------------------------
subroutine phy_growth_rate(self,Q_N,Q_P,F_T,F_co2,F_par,P_growth_rate)
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: Q_N ! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
  real(rk),dimension(self%phyto_num),intent(in) :: Q_P ! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
  real(rk),dimension(self%phyto_num),intent(in) :: F_T ! Temperature dependency for phytoplankton
  real(rk),dimension(self%phyto_num),intent(in) :: F_co2 ! CO2 forcing
  real(rk),dimension(self%phyto_num),intent(in) :: F_par ! PAR forcing
  real(rk),dimension(self%phyto_num),intent(out) :: P_growth_rate ! Phytoplankton growth rate, d^-1

  real(rk),dimension(self%phyto_num)::f_nut,r,n,g_N,mu_max,q_Nl,q_Pl

  if (self%convert_mu .eqv. .true.) then
    mu_max = self%mumax*(self%QN_max/(self%QN_max-self%QN_min))*(self%QP_max/(self%QP_max-self%QP_min)) !OG
  else
    mu_max = self%mumax !OG
  end if

  q_Nl = (Q_N-self%QN_min)/Q_N !OG
  q_Pl = (Q_P-self%QP_min)/Q_P !OG

  select case(self%Nut_lim)
    case(1)
      r = q_Pl/q_Nl
      n = self%n_star*(1.0_rk+q_Nl)
      g_N = (r-r**(1.0_rk+n))/(1.0_rk-r**(1.0_rk+n))
      f_nut = q_Nl*g_N
    case(2)
      f_nut = min(q_Nl,q_Pl)
    case(3)
      f_nut = q_Nl*q_Pl/(q_Nl+q_Pl)
    case(4)
      f_nut = q_Nl*q_Pl
    case default
      f_nut = 1.0_rk
  end select

  P_growth_rate =  mu_max*f_nut* F_T(1)*F_co2*F_par

  return
end subroutine phy_growth_rate

!------------------------------------------------------------------------------
subroutine aggr_rate(self,Phy,Q_N,D_N,aggr)! Aggregation rate, d^-1
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: Phy ! Phytoplankton biomass concentration, mmol-C m^-3
  real(rk),dimension(self%phyto_num),intent(in) :: Q_N ! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
  real(rk),intent(in) :: D_N ! Nitrogen content of detritus concentration,  mmol-N m^-3
  real(rk),intent(out) :: aggr ! actual aggregation rate

  aggr = self%A_star_opt*(sum(Phy*Q_N)+D_N)

  return
end subroutine aggr_rate

!------------------------------------------------------------------------------
subroutine N_uptake(self,N,Q_N,F_T,par,uptake_rate_N)
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),intent(in) ::  N! Nitrogen concentration, mmol-N m^-3
  real(rk),dimension(self%phyto_num),intent(in) :: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
  real(rk),dimension(self%phyto_num),intent(in) :: F_T! Temperature dependency
  real(rk),intent(in) :: par !PAR from data
  real(rk),dimension(self%phyto_num),intent(out) :: uptake_rate_N! Phytoplankton nitrogen uptake rate, mol-N mol-C^-1 d^-1

  real(rk),dimension(self%phyto_num):: nom_N,dom_N,q
  uptake_rate_N(:) = 0.0_rk

  if (par>=0.0_rk) then  !No uptake during night
    nom_N = self%vN_max*self%N_affin*N !OG
    dom_N = self%vN_max+self%N_affin*N !OG
    q = (self%QN_max-Q_N)/(self%QN_max-self%QN_min)
    where(q .le. 0.0_rk) q = 0.0_rk
    uptake_rate_N = (nom_N/dom_N)*sqrt(F_T(1))*q
  end if

  return
end subroutine N_uptake

!------------------------------------------------------------------------------
subroutine P_uptake(self,P, Q_P,F_T,par,uptake_rate_P)
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),intent(in) ::  P! Phsphorous concentration, mmol-P m^-3
  real(rk),dimension(self%phyto_num),intent(in) :: Q_P! Phytoplankton intracellular phosphorous cell quota, mol-P mol-C^-1
  real(rk),dimension(self%phyto_num),intent(in) :: F_T! Temperature dependency for phytoplankton
  real(rk),intent(in) :: par !PAR from data
  real(rk),dimension(self%phyto_num),intent(out) :: uptake_rate_P! Phytoplankton nutrient uptake rate, mol-P mol-C^-1 d^-1

  real(rk),dimension(self%phyto_num) :: nom_P,dom_P,q
  uptake_rate_P(:) = 0.0_rk

  if (par>=0.0_rk) then !No uptake during night
      nom_P = self%vP_max*self%P_affin*P !OG
      dom_P = self%vP_max+self%P_affin*P !OG
      uptake_rate_P = (nom_P/dom_P)*sqrt(F_T(1))
  end if

  return
end subroutine P_uptake

!------------------------------------------------------------------------------
subroutine sink_rate(self,Q_N,Q_P,mixl,sinking)
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: Q_N  ! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
  real(rk),dimension(self%phyto_num),intent(in) :: Q_P  ! Phytoplankton intracellular phosphorous cell quota, mol-P mol-C^-1
  real(rk),dimension(self%phyto_num),intent(out) :: sinking ! sinking rate, d^-1

  real(rk),dimension(self%phyto_num) :: physiol,qN,qP
  real(rk),intent(in) :: mixl ! mixing layer depth
  sinking = 0.0_rk

  qN = (Q_N-self%QN_min)/(self%QN_max-self%QN_min) !OG
  qP = (Q_P-self%QP_min)/(self%QP_max-self%QP_min) !OG
  physiol = exp(-0.5*((qP*qN)*16.0)**2) !Todo: Why?
  sinking = physiol*exp(0.5_rk*self%log_ESD)* 0.06_rk/mixl

  return
end subroutine sink_rate

!------------------------------------------------------------------------------
subroutine respiration(self,N_uptake,R_N)
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: N_uptake ! Nitrogen uptake rate, mol-N mol-C^-1 d^-1
  real(rk),dimension(self%phyto_num),intent(out) :: R_N ! Phytoplankton respiration rate, d^-1

  R_N = N_uptake*self%mol_ratio

  return
end subroutine respiration

!------------------------------------------------------------------------------
subroutine dD_N_dt(self, Phy, Q_N, D_N,aggregation,F_T,mixl,grazing_forc,dD_N_dt_out)
  implicit none
  !Detritus concentration over time, mmol-N m^-3 d^-1
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
  real(rk),dimension(self%phyto_num),intent(in) :: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^-1
  real(rk),intent(in) ::  D_N! Nitrogen content of detritus concentration,  mmol-N m^-3
  real(rk),intent(in) ::  aggregation! aggregation rate, d^-1
  real(rk),dimension(2),intent(in) :: F_T! Temperature dependency
  real(rk),intent(in) :: mixl
  real(rk),dimension(self%phyto_num),intent(in) :: grazing_forc! Grazing forcing, mmol-C m^-3 d^-1
  real(rk),intent(out) :: dD_N_dt_out

  real(rk) :: mtotal
  real(rk),dimension(self%phyto_num) :: source,gr,loss

  mtotal = 0.0_rk
  gr = (1.0_rk-self%y)*grazing_forc*Q_N
  loss = (self%frac_md*self%m + aggregation)*Phy*Q_N
  dD_N_dt_out = sum(loss+gr)-(self%r_dn*F_T(1)+(self%det_sink_r/mixl))*D_N

  return
end subroutine dD_N_dt
!------------------------------------------------------------------------------

subroutine dD_P_dt(self,Phy,Q_P, D_P,aggregation,F_T,mixl,grazing_forc,dD_P_dt_out) !Detritus concentration over time, mmol-P m^-3 d^-1
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
  real(rk),dimension(self%phyto_num),intent(in) :: Q_P! Phytoplankton intracellular phosphorous cell quota, mol-P mol-C^-1
  real(rk),intent(in) :: D_P! Phosphorous content of detritus concentration,  mmol-P m^-3
  real(rk),intent(in) :: aggregation! aggregation rate, d^-1
  real(rk),dimension(2),intent(in) :: F_T! Temperature dependency
  real(rk),dimension(self%phyto_num),intent(in) :: grazing_forc! Grazing forcing, mmol-C m^-3 d^-1
  real(rk),intent(in) :: mixl !mixed layer depth
  real(rk),intent(out) :: dD_P_dt_out

  real(rk) :: mtotal
  real(rk),dimension(self%phyto_num) :: source,gr,loss

  mtotal = 0.0_rk
  gr = (1.0_rk-self%y)*grazing_forc*Q_P
  loss = (self%frac_md*self%m + aggregation)*Phy*Q_P
  dD_P_dt_out=sum(loss+gr)-(self%r_dn*F_T(1)+(self%det_sink_r/mixl))*D_P
  return

end subroutine dD_P_dt

!------------------------------------------------------------------------------
subroutine dN_dt(self,N_uptake, Phy, D_N,grazing_forc,Q_N,F_T,N,dN_dt_out) ! Change of nitrogen concentration over time, mmol-N m^-3 d^-1
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: N_uptake! Phytoplankton nitorgen uptake rate, mol-N mol-C^-1 d^-1
  real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
  real(rk),intent(in) ::  D_N! Nitrogen content of detritus concentration,  mmol-N m^-3
  real(rk),dimension(self%phyto_num),intent(in) :: grazing_forc! Grazing forcing, mmol-C m^-3 d^-1
  real(rk),dimension(self%phyto_num),intent(in) :: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^1
  real(rk),dimension(2),intent(in) :: F_T! Temperature dependency
  real(rk),intent(in) ::  N! Nitrogen concentration, mmol-N m^-3
  real(rk),intent(out) ::  dN_dt_out! Nitrogen concentration, mmol-N m^-3

  real(rk),dimension(self%phyto_num) :: up

  up = N_uptake*Phy
  dN_dt_out = self%r_dn*D_N*F_T(1)-sum(up)!+sum(gr(:))

  return
end subroutine dN_dt

!------------------------------------------------------------------------------
subroutine dP_dt(self,P_uptake, Phy, D_P,grazing_forc,Q_P,F_T,P,dP_dt_out) ! Change of nutrient concentration over time
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: P_uptake! Phytoplankton phosphorous uptake rate, mol-P mol-C^-1 d^-1
  real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
  real(rk),intent(in) :: D_P! Phosphorous content of detritus concentration,  mmol-P m^-3
  real(rk),dimension(self%phyto_num),intent(in) :: grazing_forc! Grazing forcing, mmol-C m^-3 d^-1
  real(rk),dimension(self%phyto_num),intent(in) :: Q_P! Phytoplankton intracellular phosphorous cell quota, mol-P mol-C^-1
  real(rk),dimension(2),intent(in) :: F_T! Temperature dependency
  real(rk),intent(in) :: P! Phosphorous concentration, mmol-P m^-3
  real(rk),dimension(self%phyto_num) :: up
  real(rk),intent(out) ::  dP_dt_out! Nitrogen concentration, mmol-N m^-3

  up = P_uptake*Phy
  dP_dt_out = self%r_dn*D_P*F_T(1)-sum(up)

  return
end subroutine dP_dt

!------------------------------------------------------------------------------
subroutine Rel_growth_rate_sr(self,Phy,aggregation,growth_rate,grazing_forc,respiration,sinking,rel_growth_rate)
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
  real(rk),dimension(self%phyto_num),intent(out) :: rel_growth_rate
  real(rk),intent(in) ::  aggregation! Aggregation rate, d^-1
  real(rk),dimension(self%phyto_num),intent(in) :: growth_rate! Phytoplankton growth rate, d^-1
  real(rk),dimension(self%phyto_num),intent(in) :: grazing_forc! Grazing forcing, mmol-C m^-3 d^-1
  real(rk),dimension(self%phyto_num),intent(in) :: respiration! Phytoplankton respiration rate, d^-1
  real(rk),dimension(self%phyto_num),intent(in) :: sinking! Phytoplankton sinking rate, d^-1

  rel_growth_rate = growth_rate - respiration - sinking - self%m - aggregation - grazing_forc/(eps+Phy)

  return
end subroutine Rel_growth_rate_sr

!------------------------------------------------------------------------------
real(rk) function chl_a(self,Phy,Q_N,par)
implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),intent(in) :: par
  real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
  real(rk),dimension(self%phyto_num),intent(in) :: Q_N! Phytoplankton intracellular nitrogen cell quota, mol-N mol-C^1

  real(rk) :: phyto_conc_tot

  phyto_conc_tot = sum(Phy*Q_N) !* (1._rk-par/self%I_opt) !OG
  chl_a=phyto_conc_tot*self%chla_to_T_PhyN

return
end function chl_a
!-----------------------------------------------------------------
subroutine mean_cell_size(self,Phy,mean)
implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num), intent(in) ::Phy
  real(rk),dimension(self%phyto_num) :: mean_nom
  real(rk) :: Phy_tot
  real(rk), dimension(self%phyto_num) :: Phy_tmp
  integer :: i
  real(rk),intent(out) :: mean
!        :param Phy: Phytoplankton biomass concentration, mmol-C m^-3
!        :return: community mean cell size, log_e ESD (mu m)
!        """
mean_nom=0.0_rk
mean=0.0_rk
Phy_tmp=Phy
	do i=1,self%phyto_num
	!Larger Diatoms (L>3.5) were counted seperately
	        if(self%log_ESD(i)<self%log_ESD_crit) then! .and. self%log_ESD(i)/= 2.5) then
	        !Error in data for nan IV (2.6)
		  mean_nom(i)=Phy(i)*exp(self%log_ESD(i))
		!else if (self%log_ESD(i)==self%log_ESD_crit) then
		!  mean_nom(i)=0.5*Phy(i)*self%log_ESD(i)
		!  Phy_tmp(i)=0.5*Phy_tmp(i)
		else
		  Phy_tmp(i)=0.0_rk
		end if
	end do
        if (sum(Phy_tmp) /= 0.0_rk) mean = log(sum(mean_nom)/sum(Phy_tmp))
        return
end subroutine mean_cell_size

!------------------------------------------------------------------------------

real(rk) function size_diversity(self,Phy,mean)
  implicit none
  class (type_hzg_mspec),intent(in) :: self
  real(rk),dimension(self%phyto_num),intent(in) :: Phy! Phytoplankton biomass concentration, mmol-C m^-3
  ! :return: community size diversity, log_e ESD (mu m)^2
  integer::i
  real(rk),dimension(self%phyto_num) :: div_nom
  real(rk) :: mean
  !allocate(div_nom(self%phyto_num))
  do i=1,self%phyto_num
    div_nom(i)=(self%log_ESD(i)-mean)**2*Phy(i)
  end do
  size_diversity=sum(div_nom)/sum(Phy)

  return
end function size_diversity

!------------------------------------------------------------------------------
real(rk) function allometries_esd(beta,alpha,s)
  real(rk), intent(in) :: beta! Interception for trait
  real(rk), intent(in) :: alpha! Size scaling exponent for trait
  real(rk), intent(in) :: s! Equivalent spherical diamater, mu m
  ! :return: trait=exp^(log_e(beta)+alpha*log_e(ESD))
  allometries_esd =beta*exp(alpha*s)

  return
end function allometries_esd

!------------------------------------------------------------------------------

subroutine bgc_parameters(self,s, bgc_params)
  implicit none!    :return: Phytoplankton eco-physiological traits
  class (type_hzg_mspec),intent(in) :: self
  real(rk),intent (in) :: s! Equivalent spherical diamater, mu m
  real(rk),dimension(11),intent (out)::bgc_params
  real(rk) :: mu_max, Qmin_N, Qmax_N, vmax_N, Kn_N, affinity_N
  real(rk) :: Qmin_P, Qmax_P, vmax_P, Kn_P, affinity_P,tmp

  !! Nonlinear mumax
  mu_max=self%pars(3)*min(self%pars(1)*exp(s*self%pars(2)), self%pars(4)*exp(s*self%pars(5)))

 !   tmp=self%pars(1)*exp(2.29_rk*self%pars(2))
 !   mu_max=self%pars(3)*min(self%pars(1)*exp(s*self%pars(2)),tmp*exp((s-tmp)*self%pars(5)))

  !linear mumax
  !mu_max=allometries_esd(self%pars(26),self%pars(27),s)

  Qmin_N=allometries_esd(self%pars(6),self%pars(7),s)
  Qmax_N=allometries_esd(self%pars(8),self%pars(9),s)
  vmax_N =allometries_esd(self%pars(10),self%pars(11),s)

  !if (log(s)<2.3_rk) then
  !    affinity_N=self%pars(12)*dexp(self%pars(13)*2.3_rk)
  !else

  !end if
  !Kn_N =allometries_esd(self%pars(14),self%pars(15),s)
  affinity_N=allometries_esd(self%pars(12),self%pars(13),s)
  Qmin_P=allometries_esd(self%pars(16),self%pars(17),s)
  Qmax_P=allometries_esd(self%pars(18),self%pars(19),s)
  vmax_P =allometries_esd(self%pars(20),self%pars(21),s)
  !Kn_P =allometries_esd(self%pars(22),self%pars(23),s)
  affinity_P=allometries_esd(self%pars(24),self%pars(25),s)

  bgc_params(1) = mu_max
  bgc_params(2) = Qmin_N
  bgc_params(3) = Qmax_N
  bgc_params(4) = vmax_N
  bgc_params(5) = 0.0!Kn_N
  bgc_params(6) = affinity_N

  bgc_params(7) = Qmin_P
  bgc_params(8) = Qmax_P
  bgc_params(9) = vmax_P
  bgc_params(10) = 0.0!Kn_P
  bgc_params(11) = affinity_P

  return
end subroutine bgc_parameters

!------------------------------------------------------------------------------
subroutine convert_BGCparams(beta,alpha,a_carbon,b_carbon,ESD_beta,ESD_alpha)
  implicit none
  real(rk),intent(in) :: beta ! :param beta_params: Interception for trait
  real(rk), intent(in) :: alpha ! :param alpha_params: Size scaling exponent for trait
  real(rk), intent(in) :: a_carbon! size scaling exponent for cell carbon content
  real(rk), intent(in) :: b_carbon! Interception for cell carbon content
  real(rk), intent(out) :: ESD_beta,ESD_alpha

  real(rk) :: beta_cell_moleC

  beta_cell_moleC=(beta/b_carbon)*12._rk*10_rk**6!  #moleX/moleC
  ESD_beta=beta_cell_moleC*pi_over_six**(alpha-a_carbon)
  ESD_alpha=3.0_rk*(alpha-a_carbon)

  return
end subroutine convert_BGCparams

!------------------------------------------------------------------------------

!Umwandlsung von integer in Character Typ-----------------------------------
character*3 function int2char(i)
  implicit none
  integer, intent(in) :: i

  write(int2char,'(i3)') i
  int2char = adjustl(int2char)

end Function int2char

!------------------------------------------------------------------------------
subroutine make_phyto_allometries(self)
  !! Calculate and store all allometries for phytoplankton.
  !! This routine makes the calculations more efficient than the original code.
  !! OG 01.05.2019

  implicit none
  class (type_hzg_mspec), intent(inout) :: self
  real(rk), dimension(11) :: bgc_params
  integer :: i

  do i=1,self%phyto_num
    call bgc_parameters(self,self%log_ESD(i),bgc_params)
    self%mumax(i) = bgc_params(1)
    self%N_affin(i) = bgc_params(6)
    self%P_affin(i) = bgc_params(11)
    self%QN_min(i) = bgc_params(2)
    self%QN_max(i) = bgc_params(3)
    self%QP_min(i) = bgc_params(7)
    self%QP_max(i) = bgc_params(8)
    self%vN_max(i) = bgc_params(4)
    self%vP_max(i) = bgc_params(9)
  end do

  return
end subroutine make_phyto_allometries
!-------------
