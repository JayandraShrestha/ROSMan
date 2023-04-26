      Subroutine ResSman
      use parm
      use resapp
      IMPLICIT NONE
      INTEGER :: Cdate, DAYS
      real :: restime, SedIdx, settle_mass, settle_vol,Stlr_Vol_to_Dist,
     &Extra_Sed_Storage, Increm_Stlr_Vol,  Increm_Extra_Stlr_Vol,  
     & elevrule, volrule,Increm_Extra_Sed_Vol, newvol, devol, qdevol
      real :: vol, sed, vvr, targ, xx, flw, areah,yy1, yy2,yy3, yy4, yy5
      real :: min_req_release, act_max_storage, ressedo_flow, v4, v1, v2
	real :: san,sil,cla,sag,lag,gra,ndespill, v3
	real :: inised, finsed, setsed, remsetsed, flush_load, Res_len_act
      INTEGER :: jp, jk, jres
      
      
      jres = 0
      jres = inum1    
                       
    
!! store initial values
      vol = 0.
      sed = 0.
	inised = 0.
	finsed = 0.
	setsed = 0.
	remsetsed = 0.      
      vol = res_vol(jres)
      sed = res_sed(jres) 
      newvol= 0.0
      flush_load = 0.0
      Res_len_act=0.0
      act_max_storage = 0.0
      v1=0.
      v2=0.
      v3 = 0.
      v4 = 0.
      yy1 = 0.
      yy2 = 0.
      yy3 = 0.
      yy4 = 0.
      yy5 = 0.
      settle_mass = 0.
      settle_vol = 0.
      SedIdx = 0.
     
!!    Storage by particle sizes
      san = res_san(jres)
      sil = res_sil(jres)
      cla = res_cla(jres)
      sag = res_sag(jres)
      lag = res_lag(jres)
	gra = res_gra(jres) 
          !! calculate water level at the beginning of day    
      call interpolate_only(vea_curve(jres)%num_vea_data, 
     & vea_curve(jres)%voldata, vea_curve(jres)%elevdata,
     & vol, pre_elev(jres))
      
      !! calculate surface area for day    
      call interpolate_only(vea_curve(jres)%num_vea_data, 
     & vea_curve(jres)%voldata, vea_curve(jres)%areadata,
     & res_vol(jres), ressa)

      !! calculate water balance for day
      resev = 10. * evrsv(jres) * pet_day * ressa
      ressep = res_k(jres) * ressa * 240.
      respcp = sub_subp(res_sub(jres)) * ressa * 10.
      
      IF(Do_flushing(jres) ==1)then
      pre_flushing(jres) = Flushing(jres)
      do jp = 1, flush_data(jres)%fd
      call jtodate(iyr, iida, Cdate)
      If (Cdate==flush_data(jres)%flush_start(jp))then
          if (flush_data(jres)%minresfli(jp)>0.)then
                  fl_minresfli(jres) = flush_data(jres)%minresfli(jp)
          else
                  fl_minresfli(jres) = flush_data(jres)%minresfli(1)
          endif
          !! make loop to do more than one time flushing
      If (fl_minresfli(jres) < resflwi/86400.) Then !'Only attempt drawdown in this time period, and read data from matrix into variables, if known reservoir inflow is high enough.
        !'Whether or not to perform flushing
            Flushing(jres) = 1 !'reservoir should have flushing performed for the specified number of days, starting now
            Flushing_Event_Number(jres) = Flushing_Event_Number(jres)+ 1
            Cur_jp_Flush_data(jres) = jp
            Flushing_Drawdown_Start_Date(jres)=
     & flush_data(jres)%flush_start(jp)

            !'Flushing duration
            If (flush_data(jres)%flush_start(jp) > 0)Then
                Fl_Dur(jres) = flush_data(jres)%fdays(jp)
          flush_data(jres)%Fl_Dur_TEMP(Flushing_Event_Number(jres)) = 
     &           flush_data(jres)%fdays(jp)
            Else
                Fl_Dur(jres) = flush_data(jres)%fdays(1)                
      flush_data(jres)%Fl_Dur_TEMP(Flushing_Event_Number(jres)) = 
     &           flush_data(jres)%fdays(1)
            End If

            Flushing_Duration_Initial(jres) = fl_dur(jres)
            counter_22(jres) = jp
            If (flush_data(jres)%maxfl(jp) > 0) Then
                fl_maxfl(jres) = flush_data(jres)%maxfl(jp)
            Else
            !'User did not specify an appropriate value; Load the default value stored in the first row
                fl_maxfl(jres) = flush_data(jres)%maxfl(1)
            End If
        !'Minimum flushing flow
        If (flush_data(jres)%minflq(jp) > 0) Then
            fl_minflq(jres) = flush_data(jres)%minflq(jp)
        Else
        !'User did not specify an appropriate value; Load the default value stored in the first row
            fl_minflq(jres) = flush_data(jres)%minflq(1)
        End If

        !' Flushing maximum drawdown rate
        If (flush_data(jres)%mxddr(jp) > 0) then
            fl_mxddr(jres) = flush_data(jres)%mxddr(jp)
        Else
            fl_mxddr(jres) = flush_data(jres)%mxddr(1)
        End If
          !' Optional user input - flushing equation coefficient
          If (flush_data(jres)%coef_val(jp) > 0) then
              fl_coef(jres) = flush_data(jres)%coef_val(jp)
          Else
          !'User elects to have model determine flushing coefficient
              If (flush_data(jres)%coef_val(1) > 0) Then
                  fl_coef(jres) = flush_data(jres)%coef_val(1)
              Else
                  fl_coef(jres) = 0.
              End If
          End If
      !' Optional user input - flushing equation exponent
      If (flush_data(jres)%expo_val(jp)> 0) Then
          fl_expo(jres) = flush_data(jres)%expo_val(jp)
      Else
      !'User elects to have model determine flushing exponent
          If (flush_data(jres)%coef_val(1) > 0) Then
              fl_expo(jres)  = flush_data(jres)%expo_val(1)
          Else
              fl_expo(jres)  = 0
          End If
      End If

      EXIT
      Else
      !'Change flushing day in table to next day, as the inflow requirement for drawdown was not satisfied today.
      CALL Add1DAY(flush_data(jres)%flush_start(jp),1, 
     & flush_data(jres)%flush_start(jp))
      Exit
      End If
      Else
      !'The flushing variables will be reset = 0 once flushing is complete
      End If      !
      enddo
      else
      ! not exist flushing reservoirs
      ENDIF
      
   
        !' Before sluicing loop, see if sluicing completed in last loop, to see whether or not to begin post-sluicing refill routine.
      If (Post_Sluicing_1(jres) == 1) Then
          Post_Sluicing(jres) = 1
      End If
        !'Load in Sluicing Data, if sluicing is to occur in this reservoir
        If (Do_Sluicing(jres) == 1) Then
            pre_sluicing(jres) = Sluicing(jres)
            Do jp = 1, sluice_data(jres)%sd
                call jtodate(iyr, iida, Cdate)

          If (sluice_data(jres)%sluice_dur(jp) > 0) Then
          Sl_Duration(jres) = int(sluice_data(jres)%sluice_dur(jp))                   
          CALL Add1DAY(sluice_data(jres)%sluice_start(jp), 
     &     Sl_Duration(jres), Sl_end_date(jres))                    
                Else
                    !'Load default minimum inflow requirement from first SLUICING date (or zero if none specified)
              Sl_Duration(jres)= int(sluice_data(jres)%sluice_dur(1))
              sluice_data(jres)%sluice_dur(jp) = Sl_Duration(jres)                  
              CALL Add1DAY(sluice_data(jres)%sluice_start(jp), 
     &             Sl_Duration(jres), Sl_end_date(jres))                    
           End If

      If((Cdate >= sluice_data(jres)%sluice_start(jp)).and.
     & (Cdate <= Sl_end_date(jres))) Then
      If(cDATE == sluice_data(jres)%sluice_start(jp))then
      If (sluice_data(jres)%slminflq(jp) > 0.) Then
          sl_slminflq(jres) = sluice_data(jres)%slminflq(jp)
          Else
                            !'Load default minimum inflow requirement from first SLUICING date (or zero if none specified)
          sl_slminflq(jres) = sluice_data(jres)%slminflq(1)
      End If

      If ((resflwi/86400. > sl_slminflq(jres)) .and.
     &  (Flushing(jres) == 0)) Then !'Attempt drawdown in this time period, and read data from matrix into variables, if known reservoir inflow is high enough.
      Sluicing(jres) = 1 !'reservoir should be sluiced during this time period
      Cur_jp_Sluicing_data(jres) = jp
      !'Sluicing minimum inflow for refill
      If (sluice_data(jres)%slminfli(jp) > 0) Then
          sl_slminfli(jres) = sluice_data(jres)%slminfli(jp)
      Else
      !'User did not specify an appropriate value; Load the default value stored in the first row
      sl_slminfli(jres) = sluice_data(jres)%slminfli(1)
      EndIf
 !'Sluicing target drawdown elevation
      If (sluice_data(jres)%slmaxlvl(jp) > 0) Then
          sl_slmaxlvl(jres) = sluice_data(jres)%slmaxlvl(jp)
      Else
      !'User did not specify an appropriate value; Load the default value stored in the first row
      sl_slmaxlvl(jres) = sluice_data(jres)%slmaxlvl(1)
      End If
      !'Sluicing maximum drawdown rate
      If (sluice_data(jres)%slmxddr(jp) > 0) Then
      sl_slmxddr(jres) = sluice_data(jres)%slmxddr(jp)
      Else
      !'User did not specify an appropriate value; Load the default value stored in the first row
      sl_slmxddr(jres) = sluice_data(jres)%slmxddr(1)
      End If

      !'Sluicing maximum refill rate
      If (sluice_data(jres)%mxrefil(jp) > 0) Then
          sl_mxrefil(jres) = sluice_data(jres)%mxrefil(jp)
      Else
      !'User did not specify an appropriate value; Load the default value stored in the first row
      sl_mxrefil(jres) = sluice_data(jres)%mxrefil(1)
      End If

      !'Power production during sluicing?
      If (sluice_data(jres)%slpower(jp) > 0) Then
      sl_power(jres) = 1
      Else
      sl_power(jres) = 0
      End If
      Else
      !'Change sluicing start and stop days in table to next day, as the inflow requirement for drawdown was not satisfied today.  
      CALL Add1DAY(sluice_data(jres)%sluice_start(jp), 1, 
     &  sluice_data(jres)%sluice_start(jp))  
      Sluicing(jres) = 0 !'We are on the desired sluicing start date, but inflows are inadequate to begin sluicing, so today we aren't sluicing, which will ensure flushing can begin.
      Exit
      End If
      Else
      End If
      If (cDATE == Sl_end_date(jres)) Then
          If (sl_slminfli(jres) > 0.) Then
          If (resflwi/86400. > sl_slminfli(jres)) Then
      sluice_data(jres)%sluice_dur(jp) =sluice_data(jres)%sluice_dur(jp)
     &  + 1
      Else
      Post_Sluicing_1(jres) = 1
      If (Non_date_sluicing(jres) == 1) Then !'If sluicing is only triggered by inflow, not dates, then reset the next sluicing start date to tomorrow.
      CALL Add1DAY(Sl_end_date(jres),1,
     & sluice_data(jres)%sluice_start(jp))                                 
      sluice_data(jres)%sluice_dur(jp) = 1
      End If
      End If
      Else
      !'no minimum sluicing inflow value for refill was set by user
      Post_Sluicing_1(jres) = 1
      If (Non_date_sluicing(jres) == 1) Then !'If sluicing is only triggered by inflow, not dates, then reset the next sluicing start date to tomorrow.
      CALL Add1DAY(Sl_end_date(jres), 1, 
     & sluice_data(jres)%sluice_start(jp))  
      sluice_data(jres)%sluice_dur(jp) = 1                              
      End If
      End If
      Else
      End If
      !'Not the first or last sluicing day, but still a sluicing day.
      Sluicing(jres) = 1 !'reservoir should be sluiced during this time period
      Exit
      Else
      !' The currently searched row in the sluicing specifications matrix doesn't have a sluicing time frame within which the current date falls.
      Sluicing(jres) = 0
      !'Skip to next user-specified date for sluicing in the list.
      End If
      End Do
      Else
            !'There exist sluicing data for the reservoir, but the user has indicated that this function should not be performed during the simulation
      End If
      !'Set a variable value that will allow capacity of hydropower to be set to zero when sluicing is occurring if no power production is allowed during sluicing.
      If (Sluicing(jres) == 1) Then
            If (sl_power(jres) == 1) Then
                hydro_cal(jres) = 1.
            Else
                hydro_cal(jres) = 0.
            End If
      Else
            hydro_cal(jres) = 1.
      End If

        !'Establish which sediment management techniques will be given priority in the event that two techniques need to be simulated
        If ((pre_flushing(jres) == 1) .and. (Flushing(jres) == 1)) Then
            If (Sluicing(jres) == 1) Then
                !'Change sluicing start and stop days in table to next day, as sluicing must be delayed due to flushing    
      CALL Add1DAY(sluice_data(jres)%sluice_start
     & (Cur_jp_Sluicing_data(jres)), 1, sluice_data(jres)%sluice_start
     & (Cur_jp_Sluicing_data(jres))) 
                
            Else
            End If
            Sluicing(jres) = 0
        ElseIf ((pre_sluicing(jres) == 1).And.(Sluicing(jres) == 1))Then
            If (Flushing(jres) == 1) Then
                CALL Add1DAY(flush_data(jres)%flush_start
     &        (Cur_jp_Flush_data(jres)),1,flush_data(jres)%flush_start
     &        (Cur_jp_Flush_data(jres)))
            Else
            End If
            Flushing(jres) = 0
        Else
        End If
      
      If (Flushing(jres) == 1)  Then
          restime = 0.
          bruntrap(jres) = 0.
      Elseif (Sluicing(jres) == 1)then
      !'Apply Churchill method
      Res_len_act =(pre_elev(jres) - El_min(jres))/
     & ((FSL(jres)-El_min(jres))/Reslen(jres))
      If (Res_len_act > 0.) Then
          If (resflwi>0)then
                  SedIdx = ((vol/(resflwi/86400.))**2) /Res_len_act
                  If (vol>0.)then
          If ((800. * ((SedIdx / 3.28) ** (-0.2)) - 12.) > 100.) Then
                      bruntrap(jres) = 0.
          ElseIf ((800. * ((SedIdx / 3.28) **(-0.2)) - 12.) < 0.) Then
                      bruntrap(jres) = 1.
                  Else
           bruntrap(jres)=(100.-(800.*((SedIdx/3.28)**(-0.2))-12.))/100.
                  End If
                  
              else

                  bruntrap(jres) = 0.
              endif
              If(bruntrap(jres)<0.0)bruntrap(jres)=0. 
          else
          bruntrap(jres) = 0.
          endif
      else
          bruntrap(jres) = 0.
      endif    

      Else
      !! use Brune curve method
              IF (resflwi <= 0.0)then
                  restime= 5E+12
              else
                  restime=(vol/resflwi)/365. ! restime in day (=(Volume in m3/inflow volume in m3 for a day)/365days)
              Endif
              If(restime > 5E+12)then
                  restime= 5E+12
              endif
      bruntrap(jres) = (restime**br_m(jres))/(br_a(jres) + br_b(jres) *
     & restime + br_c(jres) * (restime**2.) + br_d(jres) *(restime**3.))
              if (bruntrap(jres) > 1.)bruntrap(jres) = 1.  !! set to nres
              if (bruntrap(jres) < 0.)bruntrap(jres) = 0.  !! set to nres
      End If
     
      settle_mass = bruntrap(jres) * ressedi*1000. ! convert Ton to Kg
      settle_vol = settle_mass /sed_density(jres)

      Call Interpolate_only(vea_curve(jres)%num_vea_data,
     & s_curve(jres)%elev_sed, vea_curve(jres)%voldata, 
     & FSL(jres),vol_FSL(jres))
      Call Interpolate_only(vea_curve(jres)%num_vea_data,
     & s_curve(jres)%elev_sed,vea_curve(jres)%voldata, 
     & MOL(jres),vol_MOL(jres))
      cap_act_res(jres) =vol_FSL(jres)-vol_MOL(jres)
      Cap_dead_res(jres) = vol_MOL(jres)

      ! Sediment_Storage_Capacity_Check
      if (settle_vol >= (cap_act_res(jres) + Cap_dead_res(jres)))then
      resfull(jres) = 1
      settle_vol = cap_act_res(jres) + Cap_dead_res(jres)
      settle_mass= settle_vol*sed_density(jres)
      cap_act_res(jres) = 0.
      Cap_dead_res(jres) = 0.
      resflwo = resflwi !Reservoir is full of sediment, so dam now becomes man-made waterfall with no storage (run-of-river, but with no routing).
      Elseif((cap_act_res(jres) + Cap_dead_res(jres))<=
     & storage_sum_full_res(jres))then
      resfull(jres) = 1
      resflwo = resflwi
      endif

      sed_stlr_mass(jres)= sed_stlr_mass(jres) + settle_mass       !! in kg 
      sed_stlr_vol(jres)= sed_stlr_mass(jres)/sed_density(jres)
      
      Stlr_Vol_to_Dist = 0.0
      Stlr_Vol_to_Dist = Stlr_Vol_to_Dist + settle_vol

      !Update Vol-Elev-Area curve given new sediment deposition
      Increm_Extra_Sed_Vol =0.0
      Extra_Sed_Storage = 0.
      Do jp = 2, vea_curve(jres)%num_vea_data
      vea_curve(jres)%voldata(jp)=vea_curve(jres)%voldata(jp)-
     & Stlr_Vol_to_Dist*(s_curve(jres)%elev_dep(jp)-
     & s_curve(jres)%elev_dep(jp-1))

      If ((vea_curve(jres)%voldata(jp) - vea_curve(jres)%voldata(jp-1)) 
     & > 0. )Then
      Increm_Stlr_Vol = Stlr_Vol_to_Dist * (s_curve(jres)%elev_dep(jp)-
     & s_curve(jres)%elev_dep(jp-1))
      !No additional action required, as more storage space for water is still available at this elevation
      Else
      !Additional sediment cannot be stored within the incremental storage space available at this elevation
      Increm_Extra_Stlr_Vol = Abs(vea_curve(jres)%voldata(jp) 
     &  - vea_curve(jres)%voldata(jp-1)) + jp !We add i to ensure some volume is left over
      vea_curve(jres)%voldata(jp) = vea_curve(jres)%voldata(jp-1) + jp !We do this to ensure that some storage is available at all elevations, to be realistic.
      Extra_Sed_Storage = Extra_Sed_Storage + Increm_Extra_Sed_Vol
      Increm_Stlr_Vol = Stlr_Vol_to_Dist * (s_curve(jres)%elev_dep(jp)-
     & s_curve(jres)%elev_dep(jp-1)) - Increm_Extra_Sed_Vol
      End If
      do jk = (jp + 1), vea_curve(jres)%num_vea_data
      vea_curve(jres)%voldata(jk) = vea_curve(jres)%voldata(jk) - 
     &  Increm_Stlr_Vol
      enddo
      enddo

      Increm_Extra_Sed_Vol = 0.

      !Take the extra sediment, and store it at other elevations, starting at the bottom and working upwards.
      !First, try to store the extra sediment at a higher elevation
      If (Extra_Sed_Storage > 0.) Then
      do jp = 2, vea_curve(jres)%num_vea_data
      If((vea_curve(jres)%voldata(jp) - vea_curve(jres)%voldata(jp-1)) 
     &  > jp )Then
      Increm_Extra_Sed_Vol = Extra_Sed_Storage
      vea_curve(jres)%voldata(jp) = vea_curve(jres)%voldata(jp) 
     &    - Extra_Sed_Storage
      Extra_Sed_Storage = 0.0
      !Else
      If((vea_curve(jres)%voldata(jp) - vea_curve(jres)%voldata(jp-1))
     &  < 0.0)Then
      Extra_Sed_Storage = Abs(vea_curve(jres)%voldata(jp) - 
     & vea_curve(jres)%voldata(jp-1)) + jp
      Increm_Extra_Sed_Vol= Increm_Extra_Sed_Vol - Extra_Sed_Storage
      vea_curve(jres)%voldata(jp)  = vea_curve(jres)%voldata(jp-1)  + jp
      End If
      !Else
      !'Elevation level of full of sediment, so move to next elevation
      End If

      Do jk = (jp + 1), vea_curve(jres)%num_vea_data
      vea_curve(jres)%voldata(jk) = vea_curve(jres)%voldata(jk) 
     & - Increm_Extra_Sed_Vol
      Enddo
      If (Extra_Sed_Storage == 0.0) Then
      Exit !For 'All extra sediment has been allocated to an elevation
      !Else
      End If
      Enddo
      !Else
      End If
      Stlr_Vol_to_Dist = 0.0      
      
      !! calculate concetration
    
      If (vol + resflwi+respcp - resev - ressep- settle_vol <= 0.) Then 
        res_sed(jres) = 0.0  
        
        !res_san(jres) = 1.e-6
        !res_sil(jres) = 1.e-6
        !res_cla(jres) = 1.e-6
        !res_sag(jres) = 1.e-6
        !res_lag(jres) = 1.e-6
        !res_gra(jres) = 1.e-6
     
      Else
      res_sed(jres)=(sed *vol + ressedi -settle_mass/1000.)
     & /(vol+ resflwi + respcp - resev - ressep) !! ton/m3  (settle_mass = kg)    
            
         !   res_sil(jres) = (ressili + sil * vol) 
         !& /(vol+ resflwi + respcp - resev - ressep)  
         !   res_cla(jres) = (resclai + cla * vol) 
         !& /(vol+ resflwi + respcp - resev - ressep)   
         !   res_sag(jres) = (ressagi + sag * vol) 
         !& /(vol+ resflwi + respcp - resev - ressep)    
         !   res_lag(jres) = (reslagi + lag * vol) 
         !& /(vol+ resflwi + respcp - resev - ressep)    
         !   res_gra(jres) = (resgrai + gra * vol) 
         !& /(vol+ resflwi + respcp - resev - ressep)   
         !
         !   res_sed(jres) = Max(1.e-6,res_sed(jres))
         !   res_san(jres) = Max(1.e-6,res_san(jres))
         !   res_sil(jres) = Max(1.e-6,res_sil(jres))
         !   res_cla(jres) = Max(1.e-6,res_cla(jres))
         !   res_sag(jres) = Max(1.e-6,res_sag(jres))
         !   res_lag(jres) = Max(1.e-6,res_lag(jres))
         !   res_gra(jres) = Max(1.e-6,res_gra(jres))

      End If
 
      !!Compute the new active and dead storage capacities using the E-V-A curve
      If (resfull(jres) == 0) Then
      Call Interpolate_only(vea_curve(jres)%num_vea_data,
     & s_curve(jres)%elev_sed, vea_curve(jres)%voldata, FSL(jres), 
     &  vol_FSL(jres))
      Call Interpolate_only(vea_curve(jres)%num_vea_data,
     & s_curve(jres)%elev_sed, vea_curve(jres)%voldata, MOL(jres), 
     & vol_MOL(jres))
      cap_act_res(jres) =vol_FSL(jres)-vol_MOL(jres)
      Cap_dead_res(jres) = vol_MOL(jres)

      Else
      cap_act_res(jres) = 0.
      Cap_dead_res(jres) = 0.
      End If

      !Available_Flushing_Load
      !!Before updating channel elevation, compute bottom width of trapezoidal section to be removed via flushing.
      if (Do_flushing(jres) ==1)then
      if(Flushing_Channel_Elevation(jres)<Org_Flush_Chnl_Elev(jres))then
      ! Do not remove sediment until it has accumulated up to elevation of low-level outlet
      Flushing_Removal_Volume(jres) = 0.
      else
      !sediment layer is at or above elevation of low-level outlet; flushing can proceed now
      if (Flushing_Case(jres) == 2)then
        ! All sediment can be removed
      Flushing_Removal_Volume(jres) = Flushing_Removal_Volume(jres) 
     &    + settle_vol
        else
      Daily_Sediment_layer_Depth(jres) = settle_vol / ressa_ave(jres)
      Flushing_Channel_Elevation(jres)=Flushing_Channel_Elevation(jres) 
     & + Daily_Sediment_layer_Depth(jres)
      
      Current_Flushing_Channel_Top_Width(jres)= min(W_tf(jres),
     & flush_data(jres)%W_f(1)+ 2 * (Flushing_Channel_Elevation(jres)-
     & Org_Flush_Chnl_Elev(jres))/flush_data(jres)%SS_s(1))
      
      Current_Reservoir_Sediment_Top_Width(jres) = min(W_t(jres), 
     & flush_data(jres)%W_bot(1) + 2 * (Flushing_Channel_Elevation(jres)
     &  - El_min(jres))/flush_data(jres)%SS_res(1))

        if ((cap_act_res(jres) + Cap_dead_res(jres)) / 
     &(ini_cap_act_res(jres) + ini_Cap_dead_res(jres)) >=LTCR(jres))then
        Flushing_Removal_Volume (jres) = Flushing_Removal_Volume(jres) +
     &   min(settle_vol, settle_vol*
     &   (Current_Flushing_Channel_Top_Width(jres) /
     &   Current_Reservoir_Sediment_Top_Width(jres)))
        else
            Flushing_Removal_Volume(jres) = settle_vol + 
     &       Flushing_Removal_Volume(jres)
        endif
        Sediment_Fraction_Available_Flushing(jres) =
     &   Current_Flushing_Channel_Top_Width(jres)/ 
     &   Current_Reservoir_Sediment_Top_Width(jres)

      endif
      endif
      ELSE
      ENDIF
      

      If (Flushing(jres)== 1) Then
      If (fl_dur(jres) == Flushing_Duration_Initial(jres)) Then
          If(fl_coef(jres) == 0.) Then

              Max_Flushing_Load_Removed(jres) = sed_density(jres) * 
     &           Flushing_Removal_Volume(jres) / fl_dur(jres)
              counter_23 = 1
          Else
              Max_Flushing_Load_Removed(jres)= fl_coef(jres) * 
     &         ((resflwi/86400.) ** fl_expo(jres)) * resflwi
              counter_23(jres) = 0
          EndIf

      Else
          If (counter_23 (jres)== 0) Then
              Max_Flushing_Load_Removed(jres) = fl_coef(jres) * 
     &          ((resflwi/86400.) ** fl_expo(jres)) * resflwi
          Else
          EndIf

      End If
      End If

      If (resfull(jres) == 0) Then ! if reservoir not full of sediment

          If ((Flushing(jres) == 0).and.(Sluicing(jres) == 0) .and.
     &      (Post_Sluicing(jres) == 0)) Then   
                      
              call rulecurve(jres, vea_curve(jres)%num_vea_data,
     &        vea_curve(jres)%elevdata, vea_curve(jres)%voldata,
     &        iyr, iida, elevrule, volrule)  
              
          Elseif (Flushing(jres) == 1)then
              volrule = 0. !is not from the rule curve
              elevrule = Org_Flush_Chnl_Elev(jres) !is not from the rule curve
              !Make sure that achieving elevation target does not result in exceeding the user-specified maximum drawdown rate
              If (fl_mxddr(jres) > 0.) Then
                  If ((pre_elev(jres) - elevrule) > fl_mxddr(jres)) Then
                      elevrule = pre_elev(jres) - fl_mxddr(jres)
                      ! 'Determine a storage target that corresponds to the user-specified elevation target, for times OTHER than t = 0
                  End If
              End If
      !         Call Interpolate_only(vea_curve(jres)%num_vea_data,
      !&vea_curve(jres)%elevdata,
      !& vea_curve(jres)%voldata,elevrule,volrule)!volrule is caluclated wrt to elevrule
      Elseif ((Sluicing(jres) == 1) .or. (Post_Sluicing(jres) == 1))Then

                If(Sluicing(jres) ==1)then
                    elevrule = sl_slmaxlvl(jres) !is not from the rule curve
                    !Make sure that achieving elevation target does not result in exceeding the user-specified maximum drawdown rate

               If ((pre_elev(jres) - elevrule) > sl_slmxddr(jres)) Then
                        elevrule = pre_elev(jres) - sl_slmxddr(jres)
                        ! 'Determine a storage target that corresponds to the user-specified elevation target, for times OTHER than t = 0
                    Else
                    End If
                Elseif (post_sluicing(jres)==1)then
                    ! post sluicing

              If ((elevrule - pre_elev(jres)) > sl_mxrefil(jres)) Then
                        elevrule = pre_elev(jres) + sl_mxrefil(jres)
                    else
                        original_target(jres) = elevrule
                    endif
                ENDIF
      Else
      endif
                    ! # day to the next exceeds 2 m (in either direction).            !

      if ((elevrule- pre_elev(jres)) > 2.0)then
      elevrule = pre_elev(jres)+2.0
      elseif ((pre_elev(jres)-elevrule)>2.0)then
      elevrule = pre_elev(jres)-2.0  
      endif
      
      Call Interpolate_only(vea_curve(jres)%num_vea_data,
     &vea_curve(jres)%elevdata,
     & vea_curve(jres)%voldata,elevrule,volrule)!volrule is caluclated wrt to elevrule  
      
       !! 
          newvol = vol + resflwi + respcp - resev - ressep - settle_vol
              If (newvol<= 0.0)then
                  newvol=0.0
              Endif
       call interpolate_only(vea_curve(jres)%num_vea_data,
     &        vea_curve(jres)%voldata, s_curve(jres)%elev_sed,
     &        newvol, elev(jres))   
 
             
          !!to cal total capacity of outlet at RK lvl       
          
           RK_cap(jres) = max(0., newvol)
          call interpolate_only(vea_curve(jres)%num_vea_data,
     &    vea_curve(jres)%voldata, vea_curve(jres)%elevdata,
     &    RK_cap(jres), Rk_elev(jres))
      call read_all_outlet(jres, numspill, numpool, numflush, 
     &  spillvl, Qspill, poollvl, Qpool, flushlvl, Qflush)
     
        if (Rk_elev(jres)>spillvl(1))then
        call interpolate_only(numspill, spillvl, Qspill, 
     &   Rk_elev(jres), Cap_Spill_Outlet(jres))
        else
            Cap_Spill_Outlet(jres) =0.0
        endif    
         Call Interpolate_only(vea_curve(jres)%num_vea_data,
     &vea_curve(jres)%elevdata,vea_curve(jres)%voldata,spillvl(1),yy1)
        
        Cap_Spill_Outlet(jres) = min(Cap_Spill_Outlet(jres),
     &     max(0., (RK_cap(jres) - yy1)/86400.))
                
        if (Rk_elev(jres)>poollvl(1))then
          call interpolate_only(numpool, poollvl, Qpool, 
     &    Rk_elev(jres), Cap_Hydro_Outlet(jres))          
       
        else
            Cap_Hydro_Outlet(jres)=0.0
        endif  
         Call Interpolate_only(vea_curve(jres)%num_vea_data,
     &vea_curve(jres)%elevdata,vea_curve(jres)%voldata,poollvl(1),yy2)

      Cap_Hydro_Outlet(jres) = min(Cap_Hydro_Outlet(jres), 
     &   max(0., (RK_cap(jres) - yy2)/86400.))
                
        If ((Flushing(jres) == 0).And.(Sluicing(jres) == 0)) Then !'No flushing/sluicing/venting during this time period; do not use low level gates just to achieve new storage/elevation target
                !'Do nothing; no capacity is needed for low-level gates
                Cap_low_Outlet(jres) = 0.0
      else  
                    
       call interpolate_only(numflush, flushlvl, Qflush, 
     &   Rk_elev(jres), Cap_low_Outlet(jres))   
      Endif
      
      Call Interpolate_only(vea_curve(jres)%num_vea_data,
     &vea_curve(jres)%elevdata,vea_curve(jres)%voldata,flushlvl(1),yy3)
  
        Cap_low_Outlet(jres) = min(Cap_low_Outlet(jres),
     &   max(0., Rk_cap(jres) - yy3)/86400.)  
                
           
      devol= vol+resflwi+respcp-resev-ressep-settle_vol - volrule
              If (devol<= 0.0)then
                  devol=0.0
              Endif
              qdevol=devol/86400.
       xx= 0.0       
       Call Interpolate_only(vea_curve(jres)%num_vea_data,
     & s_curve(jres)%elev_sed, vea_curve(jres)%voldata, 
     & FSL(jres),xx)
       v2 =max(0., (newvol-xx)/86400.)
       yy4 = newvol -xx 

      resflwo = Min(max(qdevol,0.),(Cap_Hydro_Outlet(jres)
     &  *hydro_cal(jres) + Cap_Spill_Outlet(jres)+Cap_low_Outlet(jres)), 
     &  Max(0., vol+resflwi+respcp-resev-ressep)/86400. )*86400.
      Call Interpolate_only(vea_curve(jres)%num_vea_data,
     & s_curve(jres)%elev_sed, vea_curve(jres)%voldata, 
     & FSL(jres),act_max_storage) 
       min_req_release = max(0.,(vol + resflwi+respcp-resev-ressep - 
     & act_max_storage)/86400.)       
       resflwo = max(resflwo, min_req_release ) 
      pre_resflwo(jres) = resflwo ! store outflow for next day
      !! Final volume
              res_vol(jres) =vol + respcp + resflwi - resev
     &        - ressep - resflwo - settle_vol  
              
       !'Begin reservoir operations routines that will distribute the computed release among the outlets
              If (Sluicing(jres) == 1) Then
                     If (Sl_power(jres) == 1) Then
                       
                         !'Sluicing will be performed, and power will be produced during this process, so water should be first released through power outlets, and then through mid-levels
                      If (resflwo/86400. <=Cap_Hydro_Outlet(jres)) Then
                             !'No spilling over the spillway needs to occur
          QDesign(jres) = Min(resflwo/86400., Cap_Hydro_Outlet(jres))
                         Else
                             QDesign(jres) = Cap_Hydro_Outlet(jres)
                             !'In a reservoir that only generates hydropower, if the release flow rate is higher than the hydropower plant capacity, this flow is probably being spilled through spillway, but need to see if water elevation is high enough for spilling
                         End If 
                      else 
          Q_low_outlet(jres) = Min(resflwo/86400.,Cap_low_Outlet(jres))
                      endif
                         
                     
            else
                If(resflwo/86400.<=Cap_Hydro_Outlet(jres))then
              Qdesign(jres) = Min(resflwo/86400.,Cap_Hydro_Outlet(jres))
                else
              Qdesign(jres) = Cap_Hydro_Outlet(jres)
                endif
                
            endif
            
            If(sed_stlr_Vol(jres) >0.75*ini_cap_act_res(jres))then ! sed vol is > 75% storage cap then cant generate hydropower
                Qdesign(jres) = 0.0
            endif
            
            
            If (Flushing(jres) == 1) Then
            Q_low_outlet(jres) = Min(Max(resflwo/86400.-  Qdesign(jres),
     &        0.), Cap_low_Outlet(jres))
                 Else
            End If                         
              
              
      else
       !! do nothing reservoir is full
          Qdesign(jres)=0.0
          Cap_low_Outlet(jres) =0.0
      ENDIF    
      
      
      If (Flushing(jres) == 1)then
          call interpolate_only(vea_curve(jres)%num_vea_data,
     &    vea_curve(jres)%voldata, vea_curve(jres)%elevdata,
     &    res_vol(jres), elev_pre_flush(jres))          
      
      ENDIF
      
          !vol check
          if (res_vol(jres)<0.0)then
          res_vol(jres)=0.0
          resflwo = max(0., (vol + respcp + resflwi - resev- ressep))
          endif

          !Compute elevation of settled sediment volume, as drawdown will not need to go below this value.
          Flushing_occurred(jres) = 0
          If (Flushing(jres) == 1) Then
          ! Call Main_Flushing_Code

          If (Drawdown_End_Date(jres) /= 0) Then

              Call xDAYS(Drawdown_End_Date(jres), Cdate, DAYS) ! caluculate number of days betn current date and drwadown end dateMax_Flushing_Load_Removed
              Time_Elapsed_Since_Drawdown(jres) = days + 1
          Else
          Call xDAYS(Flushing_Drawdown_Start_Date(jres), Cdate, DAYS)
              Time_Elapsed_During_Drawdown(jres) = days + 1
          End If


          !# Main flushing sediment mass balance routine.
      If((pre_elev(jres) <= fl_maxfl(jres) ).and.(elev_pre_flush(jres)
     &         <= fl_maxfl(jres) )) Then
          If (Drawdown_End_Date(jres) == 0) Then  ! Record that drawdown is completed so we know how long we have been drawn down, but only record once, so we don't write over value each time
              Drawdown_End_Date(jres) = Cdate
              Time_Elapsed_During_Drawdown(jres) = 0 ! Variable is just used to make sure willingness to wait for drawdown has not been exceeded, so here we eliminate that constraint since drawdown has been achieved.
          End If

      if(((resflwo/86400.-fl_minflq(jres))>= 0.).or.((fl_minflq(jres)
     &         - resflwo/86400.) <= (0.2* fl_minflq(jres))))then
          Flushing_occurred(jres) = 1
          fl_dur(jres) =  fl_dur(jres) - 1
          Flushing_Load_Removed_Daily(jres) = min(sed_stlr_mass(jres), 
     &    Max_Flushing_Load_Removed(jres))
          flushed_load = Flushing_Load_Removed_Daily(jres)
          sed_stlr_mass(jres)= sed_stlr_mass(jres) - flushed_load  
          sed_stlr_Vol(jres) = sed_stlr_mass(jres)/ sed_density(jres)
          Flushing_Start_Date(jres) = cdate

        if (fl_dur(jres) == 0)then            
        Flushing(jres) = 0
        Flushing_End_Date(jres) = Cdate
        Call xDAYS(Flushing_Drawdown_Start_Date(jres), 
     &   Flushing_Start_Date(jres), days)
        counter_1(jres) = (Flushing_Start_Date(jres) -
     &   Flushing_Drawdown_Start_Date(jres)) + 1
        counter_1(jres) = days + 1
        Call xDAYS(Flushing_Start_Date(jres), cdate, days)
        counter_1(jres) = days + 1
      flush_data(jres)%Flushing_Total_Flows(Flushing_Event_Number
     & (jres))= 0
        counter_23(jres) = 0
        Flushing_Removal_Volume(jres) = 0.
        Flushing_Start_Date(jres) = 0
        Flushing_Drawdown_Start_Date(jres) = 0
        Drawdown_End_Date(jres) = 0
        Time_Elapsed_Since_Drawdown(jres) = 0
      endif
      flush_load = Flushing_Load_Removed_Daily(jres)/1000.  ! convert ot Ton
      else
      if ((Time_Elapsed_During_Drawdown(jres) >= 
     &  Flushing_Drawdown_Willingness_To_Wait(jres)).and.
     & (Flushing_Drawdown_Willingness_To_Wait(jres) > 0).and.
     & (Time_Elapsed_During_Drawdown(jres) > 0))then
      drawdown_terminated(jres) = 1
      flushing(jres) = 0
      endif
      endif
      else
      if((Time_Elapsed_During_Drawdown(jres) >= 
     & Flushing_Drawdown_Willingness_To_Wait(jres)).and.
     & (Flushing_Drawdown_Willingness_To_Wait(jres) > 0).and.
     & (Time_Elapsed_During_Drawdown(jres) > 0))then
      drawdown_terminated(jres) = 1
      flushing(jres) = 0
      endif
      endif

      if(flushing(jres) == 1)then
      if((Time_Elapsed_Since_Drawdown(jres) >= 
     & Flushing_Willingness_To_Wait(jres)) .and.
     & (Flushing_occurred(jres) == 1))then

      else
      if((Time_Elapsed_Since_Drawdown(jres) > 0).and.
     & (Time_Elapsed_Since_Drawdown(jres) >=
     & Flushing_Willingness_To_Wait(jres)))then
      Flushing(jres) = 0
      Time_Elapsed_Since_Drawdown(jres) = 0
      Drawdown_End_Date(jres) = 0
      Flushing_End_Date(jres) = Cdate

      if(( flush_data(jres)%fl_dur_temp(Flushing_Event_Number(jres))-
     &     fl_dur(jres))>0)then

      Call xDAYS(Flushing_Drawdown_Start_Date(jres), 
     &  Flushing_Start_Date(jres), counter_1(jres))
      counter_1(jres) = counter_1(jres) + 1

      Call xDAYS(Flushing_Start_Date(jres), cdate, counter_1(jres))
      counter_1(jres) = counter_1(jres) + 1

      else
      Call xDAYS(Flushing_Drawdown_Start_Date(jres), cdate, 
     &  counter_1(jres))
      counter_1(jres) = counter_1(jres) + 1

      endif

       flush_data(jres)%Flushing_Total_Flows(Flushing_Event_Number
     &  (jres)) = 0.  !# Reset so array can be used next time for Flushing
      counter_23(jres) = 0
      Flushing_Removal_Volume(jres) = Flushing_Removal_Volume(jres) - 
     &(1./ sed_density(jres)) *(flush_data(jres)%fl_dur_temp
     & (Flushing_Event_Number(jres))- Fl_Dur(jres)) * 
     & Flushing_Load_Removed_Daily(jres)  !# Not all of the sediment in the flushing channel was removed, since flushing was cut short by operator.
      Flushing_Start_Date(jres) = 0
      Flushing_Drawdown_Start_Date(jres) = 0
      ENDIF
      End If
      else
      endif


      if ((drawdown_terminated(jres)== 1) .and.(flushing(jres)== 0))then
      Time_Elapsed_Since_Drawdown(jres) = 0
      Time_Elapsed_During_Drawdown(jres) = 0
      Drawdown_End_Date(jres) = 0
      Flushing_End_Date(jres) = 0
       flush_data(jres)%Flushing_Total_Flows(Flushing_Event_Number
     &   (jres)) = 0  !# Reset so array can be used next time for Flushing
      counter_23(jres) = 0
      Flushing_Start_Date(jres) = 0
      Flushing_Drawdown_Start_Date(jres) = 0
      ENDIF
      drawdown_terminated(jres) = 0   ! # zero out variable for use next time

      else
      EndIf
      !! End of main Flushing code
      !'  Adjust E-V-A curve for sediment removal by flushing
      If (Flushing_occurred(jres) == 1) Then
      Do jp = 2, vea_curve(jres)%num_vea_data
          
          vea_curve(jres)%voldata(jp) = vea_curve(jres)%voldata(jp) + 
     &     Flushed_Load* (1. / sed_density(jres)) * 
     &     (s_curve(jres)%elev_rem(jp)- s_curve(jres)%elev_rem(jp-1))
          do jk =jp + 1, vea_curve(jres)%num_vea_data
              vea_curve(jres)%voldata(jk) = vea_curve(jres)%voldata(jk) 
     &         + Flushed_Load * (1. / sed_density(jres)) *
     &         (s_curve(jres)%elev_rem(jp)-s_curve(jres)%elev_rem(jp-1))
          end do
      End do
      Call Interpolate_only(vea_curve(jres)%num_vea_data,
     & vea_curve(jres)%elevdata, vea_curve(jres)%voldata, FSL(jres), 
     & vol_FSL(jres))
       Call Interpolate_only(vea_curve(jres)%num_vea_data,
     & vea_curve(jres)%elevdata, vea_curve(jres)%voldata, MOL(jres), 
     &  vol_MOL(jres))
      cap_act_res(jres) =vol_FSL(jres)-vol_MOL(jres)
      Cap_dead_res(jres) = vol_MOL(jres) 

      End If
      
!! sediment balance     
      ressedo_flow = res_sed(jres)*resflwo !! sediment outflow with water outflow Ton
      ressedo =flush_load + ressedo_flow !! flushing removal + sediment outflow with water outflow in Ton
       ! Final sediment concentration
      v4 = res_sed(jres) 
      
      if (res_vol(jres)>0.)then
      res_sed(jres) = (sed*vol + ressedi- (settle_mass/1000.)
     & - (v4*resflwo))/res_vol(jres)
      else
         res_sed(jres) = 0.0
      endif
      
         !! net change in amount of sediment (Ton) in water of reservoir for day
      ressedc = vol * sed + ressedi - ressedo_flow  - settle_mass/1000.
      if(ressedc<0.)then
          ressedc = 0.0
      endif
      
          !' Compute hydropower output for period of interest
        if(incap(jres)>0.)then
            if((elev(jres)-tlvl(jres))>0.)then
      MaxQdesign(jres)=incap(jres)*1000./(9.81*eff(jres)*
     &   (elev(jres)-tlvl(jres))*(1-hl(jres)))
            else
                MaxQdesign(jres)= 0.0
            endif
      else
       MaxQdesign(jres)= 0.0
      endif
      
            
      
        power(jres) = min(incap(jres),(9.81*eff(jres)*(elev(jres)
     &   -tlvl(jres))*(1.-hl(jres))*Qdesign(jres)/1000.0))    !! MW
        energy(jres)=power(jres)*24.0/1000.0 !gwh
      

      !! update final reservoir elevation and area      
      call interpolate_only(vea_curve(jres)%num_vea_data,
     &        vea_curve(jres)%voldata, s_curve(jres)%elev_sed,
     &        res_vol(jres), elev(jres))
      
      call interpolate_only(vea_curve(jres)%num_vea_data, 
     & vea_curve(jres)%voldata, vea_curve(jres)%areadata,
     & res_vol(jres), ressa)
      
             !!' If sluicing is being conducted, check the storage or elevation target error
        If (Post_Sluicing(jres) == 1) Then
          If (volrule > 0.) Then
          If (100. * ((res_vol(jres) - volrule) / volrule) >= -5.) Then
                    Post_Sluicing(jres) = 0
                    Post_Sluicing_1(jres) = 0
                Else
                End If
            Else
            End If
        Else
        End If
      
      end subroutine ResSman     
