MODULE hgrctl
  USE kinds
  USE param
  USE common
  USE functions
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------


  IMPLICIT NONE

CONTAINS

  SUBROUTINE hgr_ctl
    !!C---------------------------------------------------------------------
    !!!
    !!!                       ROUTINE hgr_ctl
    !!!                     ******************
    !!!
    !!!  Purpose :
    !!!  ---------
    !!!	Controle des caracteristiques des fonctions servant a construire
    !!!	la grille du modele global.
    !!!
    !!
    !!   Method :
    !!   --------
    !!	La grille a construire dans l hemisphere nord est determinee a
    !!	partir des trois fonction de base: FSY0, FSA, FSB exprimees dans
    !!	le plan stereographique polaire.
    !!      Ces fonctions doivent verifier les contraintes suivantes :
    !!
    !!       -conditions aux limites:
    !!           FSA(1)       = req     L ellipse est confondue avec le
    !!                                  cercle equateur de la grille
    !!           FSB(1)       = req     equateur
    !!           FSB(JPNORD)  = 0       L ellipse est reduite a une droite
    !!           FSY0(1)      = 0       Le centre de l ellipse part du
    !!                                  pole nord geographique
    !!           FSY0(JPNORD) = rpol    jusqu au centre de la grille
    !!
    !!       -contraintes de non-recouvrement:
    !!           FSDA         < 0       FSA decroissante
    !!           FSDB         < 0       FSB decroissante
    !!           abs(FSDY0)<abs(FSDA)
    !!
    !!	Si ces contraintes ne sont pas respectees, le programme	stoppe.
    !!
    !!
    !!   Modifications:
    !!   --------------
    !!       Original  : 92-10 (G. MADEC)
    !!!---------------------------------------------------------------------
    REAL(wp) :: zzerom, zj, zy1, zy1d, zy2, zy2d, za1, zb1, za2, zb2, za, zb, zy

    INTEGER :: istop, jj, ideca, idecb, idecy

    !!!---------------------------------------------------------------------
    !!!  OPA7, LODYC (1/7/92)
    !!!---------------------------------------------------------------------
!!$    C
    zzerom = 1.E-14
    istop = 0
!!$    C
!!$    C ================================================
!!$    C impression des fonctions de base et leur derivee
!!$    C ================================================
!!$    C     
    if ( ndebug .GE. 1 ) then
       write(0,*)
       write(0,*) ' fonction de base    A         et         B '
       write(0,*) '---------------------------------------------'
       write(0,9001)
       do jj = 1, jpnord
          zj  =FLOAT(jj)
          zy1 =fsa(zj)
          zy1d=fsda(zj)
          zy2 =fsb(zj)
          zy2d=fsdb(zj)
          write(0,9011) jj, zy1, zy1d, zy2, zy2d
       end do
       write(0,*)
       write(0,*) ' fonction de base    PHIA      et       PHIB'
       write(0,*) '---------------------------------------------'
       write(0,9001)
       do jj = 1, jpnord
          zj  =FLOAT(jj)
          zy1 =fsphia(zj)
          zy1d=fsdphia(zj)
          zy2 =fsphib(zj)
          zy2d=fsdphib(zj)
          write(0,9011) jj, zy1, zy1d, zy2, zy2d
       end do

9001   format('    jj       PHIA    PHIA derivee', &
               '   PHIB     PHIB derivee' )
9011   format(I6,3X,2F10.3,3X,2F10.3)
    end if
!!$    C
!!$    C
!!$    C ======================
!!$    C conditions aux limites
!!$    C ======================
!!$    C
!!$    C
    za1 = fsa(1.0_wp)
    zb1 = fsb(1.0_wp)
    zy1 = fsy0(1.0_wp)
    za2 = fsa (rjpnth)
    zb2 = fsb (rjpnth)
    zy2 = fsy0(rjpnth)
    write (0,*)
    write (0,*) '    Conditions aux limites'
    write (0,*) '    ----------------------'
    write (0,*)
    write (0,*) 'Grand axe                FSA ( equateur  ) = ', za1
    write (0,*) '                         FSA (pole grille) = ', za2
    write (0,*) 'Petit axe                FSB ( equateur  ) = ', zb1
    write (0,*) '                         FSB (pole grille) = ', zb2
    write (0,*) 'Centre                   FSY0( equateur  ) = ', zy1
    write (0,*) '                         FSY0(pole grille) = ', zy2



!!$    C     write (0,*) 'Grand axe     FSAdeg ( equateur  ) = ', fsyphi(za1)
!!$    C     write (0,*) '              FSAded (pole grille) = ', fsyphi(za2)
!!$    C     write (0,*) 'RPHIA = ', fsyphi(0.5*ABS(fsphiy(50.)-fsphiy(108.)))
!!$    C     write (0,*) 'pole1 = ', (rpi/2. - 2.*ATAN(ABS(fsy0(rjpnth
!!$    C    $                                            )-fsa(rjpnth))))/rad
!!$    C     write (0,*) 'pole2 = ', (rpi/2. - 2.*ATAN(ABS(fsy0(rjpnth
!!$    c    $                                            )+fsa(rjpnth))))/rad
!!$    C
!!$    C
!!$    C     
    if ((abs(ZA1-REQ).GT.zzerom).OR.(abs(ZB2).GT.zzerom).OR. &
           (abs(ZB1-REQ).GT.zzerom).OR. &
           (abs(ZY1).GT.zzerom)          ) then
       write (0,*)
       write (0,*) 'ERREUR : les conditions aux limites ne sont pas', &
               ' respectee'
       write (0,*) '********'
       print*, (abs(ZA1-REQ).GT.zzerom), abs(ZA1-REQ)
       print*, (abs(ZB2).GT.zzerom), abs(ZB2)
       print*, (abs(ZB1-REQ).GT.zzerom), abs(ZB1-REQ)
       print*, (abs(ZY1).GT.zzerom), abs(ZY1) ! fsyo
       print*, "req = ",req
       ISTOP=1
    end if
!!$    C
!!$    C ==============================
!!$    C contrainte de non recouvrement
!!$    C ==============================
!!$    C
!!$    C     
    write(0,*)
    write(0,*) '    Contraintes de non recouvrement'
    write(0,*) '    -------------------------------'
    write(0,*)
    IDECA=0
    IDECB=0
    IDECY=0
    do JJ=1,JPNORD
       ZJ=float(JJ)
       ZA=FSDA (ZJ)
       ZB=FSDB (ZJ)
       ZY=FSDY0(ZJ)
       if (ZA.GE.0) then
          write(0,*) 'JJ=',JJ,'     FSDA(JJ)=',ZA
          IDECA=1
       end if
       if (ZB.GE.0) then
          write(0,*) 'JJ=',JJ,'     FSDB(JJ)=',ZB
          IDECB=1
       end if
       if (abs(ZY).GT.abs(ZA)) then
          write(0,FMT="('JJ=',i4,'  abs(FSDY0(JJ))=',e12.4,  &
                    & ' > abs(FSDA(JJ))=',e12.4)") JJ,ABS(ZY),abs(ZA)
          IDECY=1

          print*,'edge:',FSA(zj)+FSY0(zj)

!!$          c$$$            print*,"FSA:",FSA(zj-1),FSA(zj),FSA(zj+1)
!!$          c$$$            print*,"FSY0:",FSY0(zj-1),FSY0(zj),FSY0(zj+1)
!!$          c$$$
!!$          c$$$            print*,"FSDA:",FSA(zj+1)-FSA(zj-1)
!!$          c$$$            print*,"FSDY0:",FSY0(zj+1)-FSY0(zj-1)

       end if
    end do
    if (IDECA.EQ.0) then
       write(0,*) 'La fonction FSA est decroissante           :',&
               ' FSDA < 0'
    else
       ISTOP=1
       write(0,*) 'ERREUR : La fonction FSA n est pas decroissante'
       write(0,*) '********'
    end if
    if (IDECB.EQ.0) then
       write(0,*) 'La fonction FSB est decroissante           :',&
               ' FSDB < 0'
    else
       ISTOP=1
       write(0,*) 'ERREUR : La fonction FSB n est pas decroissante'
       write(0,*) '********'
    end if
    if (IDECY.EQ.0) then
       write(0,*) 'Il n y a pas recouvrement sur l axe Oy     :',&
               ' |FSDY0| < |FSDA|'
    else
       ISTOP=1
       write(0,*) 'ERREUR : |FSDY0| > |FSDA|'
       write(0,*) '********'
    end if
!!$    C
!!$    C ================================
!!$    C arret du programme si necessaire
!!$    C ================================
!!$    C 
!!$    C
    if(ISTOP.EQ.1) then
       write(0,*)
       write(0,*) '***********************************************'
       write(0,*) '*                                             *'
       write(0,*) '*    HGR_CTL: Pb dans les fonctions de base    *'
       write(0,*) '*                                             *'
       write(0,*) '***********************************************'
       write(0,*)
       stop
    else
       write(0,*)
       write(0,*) '    HGR_CTL: Fonctions de base controlees'
       write(0,*) '    ------------------------------------'
       write(0,*)
    endif
!!$    C
!!$    C
  END SUBROUTINE hgr_ctl

END MODULE hgrctl
